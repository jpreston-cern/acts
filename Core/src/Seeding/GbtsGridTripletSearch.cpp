// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsGridTripletSearch.hpp"

#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <utility>

namespace Acts::Experimental {

GbtsGridTripletSearch::GbtsGridTripletSearch(
    const Config& config, std::unique_ptr<const Logger> logger)
    : m_cfg(config),
      m_logger(std::move(logger)),
      m_filterLogger(m_logger->cloneWithSuffix("Filter")),
      m_seeder(m_logger->cloneWithSuffix("Finder")) {
  const std::size_t nZBins =
      m_cfg.grid.zBinEdges.empty() ? 0 : m_cfg.grid.zBinEdges.size() - 1;
  if (!m_cfg.middleRRanges.empty() && m_cfg.middleRRanges.size() != nZBins) {
    throw std::invalid_argument(
        "GbtsGridTripletSearch: middleRRanges needs one entry per z bin");
  }
}

void GbtsGridTripletSearch::findTriplets(
    const Inputs& inputs, const float bFieldInZ,
    std::vector<std::array<std::uint32_t, 3>>& triplets,
    std::vector<float>& quality, std::uint32_t* belowQuality) const {
  triplets.clear();
  quality.clear();

  CylindricalSpacePointGrid::Config gridCfg = m_cfg.grid;
  gridCfg.bFieldInZ = bFieldInZ;
  CylindricalSpacePointGrid grid(gridCfg, logger().cloneWithSuffix("Grid"));

  const auto nInputs = static_cast<std::uint32_t>(inputs.nodes.size());
  const auto stripOf = [&](const std::uint32_t node)
      -> const OuterStripSpacePointCalibrationDetails* {
    if (inputs.stripIndex.empty() || inputs.stripIndex[node] == kNoStrip) {
      return nullptr;
    }
    return &inputs.strips[inputs.stripIndex[node]];
  };

  // the position of each input relative to the offset; the grid holds the
  // input index, which the seeder's container keeps as copiedFromIndex
  std::vector<std::array<float, 2>> xy(nInputs);
  std::vector<float> radius(nInputs);
  std::vector<float> phi;
  std::vector<float> asinD0OverR;
  const float dPhiD0 = m_cfg.doubletDPhiD0Max;
  if (m_cfg.doubletDPhiCut) {
    phi.resize(nInputs);
    asinD0OverR.resize(nInputs);
  }
  for (std::uint32_t i = 0; i < nInputs; ++i) {
    const std::array<float, 4>& p = inputs.positions[inputs.nodes[i]];
    const float x = p[0] - inputs.offset[0];
    const float y = p[1] - inputs.offset[1];
    const float r = std::hypot(x, y);
    const float phiI = std::atan2(y, x);
    grid.insert(i, phiI, p[2], r);
    xy[i] = {x, y};
    radius[i] = r;
    if (m_cfg.doubletDPhiCut) {
      phi[i] = phiI;
      asinD0OverR[i] = std::asin(std::min(1.f, dPhiD0 / std::max(r, 1.f)));
    }
  }

  // with the layer connections, each bin is split into one range per layer,
  // sorted by (layer, radius) in one go; otherwise by radius
  const bool layerSearch = !inputs.allowedPairs.empty();
  const auto layerOf = [&](const SpacePointIndex i) {
    return inputs.layers[inputs.nodes[i]];
  };
  for (std::size_t i = 0; i < grid.numberOfBins(); ++i) {
    if (layerSearch) {
      std::ranges::sort(
          grid.at(i), [&](const SpacePointIndex a, const SpacePointIndex b) {
            const GbtsLayerIndex layerA = layerOf(a);
            const GbtsLayerIndex layerB = layerOf(b);
            return layerA != layerB ? layerA < layerB : radius[a] < radius[b];
          });
    } else {
      std::ranges::sort(grid.at(i),
                        [&](const SpacePointIndex a, const SpacePointIndex b) {
                          return radius[a] < radius[b];
                        });
    }
  }

  SpacePointContainer spacePoints;
  spacePoints.createColumns(
      SpacePointColumns::CopiedFromIndex | SpacePointColumns::PackedXY |
      SpacePointColumns::PackedZR | SpacePointColumns::VarianceZ |
      SpacePointColumns::VarianceR);
  if (m_cfg.useStripDetails) {
    spacePoints.createColumns(SpacePointColumns::StripCalibrationDetails);
  }
  spacePoints.reserve(grid.numberOfSpacePoints());
  std::vector<SpacePointIndexRange> binRanges;
  binRanges.reserve(grid.numberOfBins());
  // the layer ranges of bin b are layerRanges[layerRangesBegin[b]] up to
  // layerRanges[layerRangesBegin[b + 1]], in increasing layer
  struct LayerRange {
    GbtsLayerIndex layer;
    SpacePointIndexRange range;
  };
  std::vector<LayerRange> layerRanges;
  std::vector<std::uint32_t> layerRangesBegin;
  if (layerSearch) {
    layerRangesBegin.reserve(grid.numberOfBins() + 1);
  }
  for (std::size_t bin = 0; bin < grid.numberOfBins(); ++bin) {
    const auto begin = static_cast<std::uint32_t>(spacePoints.size());
    if (layerSearch) {
      layerRangesBegin.push_back(
          static_cast<std::uint32_t>(layerRanges.size()));
    }
    for (const SpacePointIndex i : grid.at(bin)) {
      if (layerSearch) {
        const auto index = static_cast<std::uint32_t>(spacePoints.size());
        if (layerRanges.size() == layerRangesBegin.back() ||
            layerRanges.back().layer != layerOf(i)) {
          layerRanges.push_back({layerOf(i), {index, index}});
        }
        layerRanges.back().range.second = index + 1;
      }
      const std::uint32_t node = inputs.nodes[i];
      const std::array<float, 4>& p = inputs.positions[node];
      const std::array<float, 2>& variance = inputs.variances[node];

      auto sp = spacePoints.createSpacePoint();
      sp.copiedFromIndex() = i;
      sp.xy() = xy[i];
      sp.zr() = std::array<float, 2>{p[2], radius[i]};
      sp.varianceZ() = m_cfg.varianceScale * variance[1];
      sp.varianceR() = m_cfg.varianceScale * variance[0];
      if (m_cfg.useStripDetails) {
        if (const OuterStripSpacePointCalibrationDetails* strip = stripOf(node);
            strip != nullptr) {
          OuterStripSpacePointCalibrationDetails details = *strip;
          details.outerCenter[0] -= inputs.offset[0];
          details.outerCenter[1] -= inputs.offset[1];
          sp.outerStripCalibrationDetails() = details;
        }
      }
    }
    binRanges.emplace_back(begin,
                           static_cast<std::uint32_t>(spacePoints.size()));
  }
  if (layerSearch) {
    layerRangesBegin.push_back(static_cast<std::uint32_t>(layerRanges.size()));
  }
  const auto layerRangesOf = [&](const std::size_t bin) {
    return std::span(layerRanges)
        .subspan(layerRangesBegin[bin],
                 layerRangesBegin[bin + 1] - layerRangesBegin[bin]);
  };
  // whether a doublet may go from the outer to the inner layer
  const auto connects = [&](const GbtsLayerIndex outer,
                            const GbtsLayerIndex inner) {
    const std::size_t pair =
        static_cast<std::size_t>(outer) * inputs.nLayers + inner;
    return inputs.allowedPairs[pair] != 0;
  };
  // the strip calibration details derived once per space point, not for every
  // doublet candidate
  if (m_cfg.useStripDetails) {
    DoubletSeedFinder::addDerivedStripCalibrationColumn(spacePoints);
  }

  // the radial range of the space points, from the bins sorted by radius
  float rMin = std::numeric_limits<float>::max();
  float rMax = std::numeric_limits<float>::lowest();
  for (const SpacePointIndexRange& range : binRanges) {
    if (range.first == range.second) {
      continue;
    }
    rMin = std::min(spacePoints[range.first].zr()[1], rMin);
    rMax = std::max(spacePoints[range.second - 1].zr()[1], rMax);
  }
  const Range1D<float> rMiddleRange(
      std::floor(rMin / 2) * 2 + m_cfg.middleRMinMargin,
      std::floor(rMax / 2) * 2 - m_cfg.middleRMaxMargin);

  // the azimuthal-swing doublet cut, on the inputs' cached azimuth and
  // asin(d0 / r) (asin(d0 / r) falls with r, so the swing between the two radii
  // is the magnitude of the difference)
  auto doubletSelection = [&](const ConstSpacePointProxy& middle,
                              const ConstSpacePointProxy& other,
                              float /*cotTheta*/, bool /*isBottomCandidate*/) {
    const float rM = middle.zr()[1];
    const float rO = other.zr()[1];
    const float rInner = std::min(rM, rO);
    const float rOuter = std::max(rM, rO);
    const auto iM = middle.copiedFromIndex();
    const auto iO = other.copiedFromIndex();
    float dPhi = phi[iO] - phi[iM];
    const float swing = std::abs(asinD0OverR[iO] - asinD0OverR[iM]);
    constexpr float pi = std::numbers::pi_v<float>;
    if (dPhi > pi) {
      dPhi -= 2.f * pi;
    } else if (dPhi < -pi) {
      dPhi += 2.f * pi;
    }
    const float bound = m_cfg.doubletDPhiConst +
                        m_cfg.doubletDPhiSlope * (rOuter - rInner) +
                        std::min(m_cfg.doubletDPhiCap, swing);
    return std::abs(dPhi) <= bound;
  };

  DoubletSeedFinder::Config bottomCfg = m_cfg.bottomDoublets;
  DoubletSeedFinder::Config topCfg = m_cfg.topDoublets;
  if (m_cfg.doubletDPhiCut) {
    bottomCfg.experimentCuts.connect(doubletSelection);
    topCfg.experimentCuts.connect(doubletSelection);
  }
  const auto bottomFinder = DoubletSeedFinder::create(
      DoubletSeedFinder::DerivedConfig(bottomCfg, bFieldInZ));
  const auto topFinder = DoubletSeedFinder::create(
      DoubletSeedFinder::DerivedConfig(topCfg, bFieldInZ));
  const auto tripletFinder = TripletSeedFinder::create(
      TripletSeedFinder::DerivedConfig(m_cfg.triplets, bFieldInZ));

  BroadTripletSeedFilter::State filterState;
  BroadTripletSeedFilter::Cache filterCache;
  TripletSeeder::Cache cache;
  const BroadTripletSeedFilter filter(m_cfg.filter, filterState, filterCache,
                                      *m_filterLogger);

  std::vector<SpacePointContainer::ConstRange> bottomRanges;
  std::vector<SpacePointContainer::ConstRange> topRanges;
  SeedContainer seeds;

  const std::vector<float>& zBinEdges = m_cfg.grid.zBinEdges;
  for (const auto [bottom, middle, top] : grid.binnedGroup()) {
    if (middle >= binRanges.size()) {
      throw std::runtime_error(
          "GbtsGridTripletSearch: the grid returned an unreasonable middle "
          "bin");
    }
    // all the middle space points of a group are in one z bin
    std::pair<float, float> middleR{rMiddleRange.min(), rMiddleRange.max()};
    if (!m_cfg.middleRRanges.empty() &&
        binRanges[middle].first != binRanges[middle].second) {
      const float zM = spacePoints.at(binRanges[middle].first).zr()[0];
      const auto edge = std::ranges::lower_bound(zBinEdges, zM);
      auto zBin = static_cast<std::size_t>(edge - zBinEdges.begin());
      zBin = zBin == 0 ? zBin : zBin - 1;
      middleR = {m_cfg.middleRRanges[zBin][0], m_cfg.middleRRanges[zBin][1]};
    }

    if (layerSearch) {
      // each middle layer with the bottom and top layers it connects to
      for (const LayerRange& m : layerRangesOf(middle)) {
        bottomRanges.clear();
        topRanges.clear();
        for (const std::size_t b : bottom) {
          for (const LayerRange& l : layerRangesOf(b)) {
            if (connects(m.layer, l.layer)) {
              bottomRanges.push_back(spacePoints.range(l.range).asConst());
            }
          }
        }
        if (bottomRanges.empty()) {
          continue;
        }
        for (const std::size_t t : top) {
          for (const LayerRange& l : layerRangesOf(t)) {
            if (connects(l.layer, m.layer)) {
              topRanges.push_back(spacePoints.range(l.range).asConst());
            }
          }
        }
        if (topRanges.empty()) {
          continue;
        }
        m_seeder.createSeedsFromGroups(
            cache, *bottomFinder, *topFinder, *tripletFinder, filter,
            spacePoints, bottomRanges, spacePoints.range(m.range).asConst(),
            topRanges, middleR, seeds);
      }
      continue;
    }

    bottomRanges.clear();
    topRanges.clear();
    std::ranges::transform(bottom, std::back_inserter(bottomRanges),
                           [&](const std::size_t b) {
                             return spacePoints.range(binRanges[b]).asConst();
                           });
    const SpacePointContainer::ConstRange middleRange =
        spacePoints.range(binRanges[middle]).asConst();
    std::ranges::transform(top, std::back_inserter(topRanges),
                           [&](const std::size_t t) {
                             return spacePoints.range(binRanges[t]).asConst();
                           });

    m_seeder.createSeedsFromGroups(
        cache, *bottomFinder, *topFinder, *tripletFinder, filter, spacePoints,
        bottomRanges, middleRange, topRanges, middleR, seeds);
  }

  triplets.reserve(seeds.size());
  quality.reserve(seeds.size());
  const bool forwardQualityCut =
      m_cfg.minQualityForward > std::numeric_limits<float>::lowest();
  std::uint32_t nBelowQuality = 0;
  for (const auto seed : seeds) {
    const auto indices = seed.spacePointIndices();
    if (m_cfg.seedQualitySelection) {
      const float q = seed.quality();
      const auto& best = filterState.bestSeedQualityMap;
      if (!(best.at(indices[0]) <= q || best.at(indices[1]) <= q ||
            best.at(indices[2]) <= q)) {
        continue;
      }
    }
    if (forwardQualityCut && seed.quality() < m_cfg.minQualityForward) {
      const auto bottom = spacePoints.at(indices[0]).zr();
      const auto top = spacePoints.at(indices[2]).zr();
      const float lever = top[1] - bottom[1];
      if (lever > 0.f && std::abs(std::asinh((top[0] - bottom[0]) / lever)) >=
                             m_cfg.minQualityForwardMinAbsEta) {
        ++nBelowQuality;
        continue;
      }
    }
    std::array<std::uint32_t, 3> triplet{};
    for (std::size_t k = 0; k < 3; ++k) {
      triplet[k] = inputs.nodes[spacePoints.at(indices[k]).copiedFromIndex()];
    }
    triplets.push_back(triplet);
    quality.push_back(seed.quality());
  }
  if (belowQuality != nullptr) {
    *belowQuality = nBelowQuality;
  }
}

}  // namespace Acts::Experimental
