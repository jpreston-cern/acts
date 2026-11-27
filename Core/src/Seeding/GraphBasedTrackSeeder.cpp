// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"

#include "Acts/Seeding/GbtsGraphBuilder.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "Acts/Seeding/detail/TripletGbtsGraphTypes.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Acts::Experimental {

GraphBasedTrackSeeder::DerivedConfig::DerivedConfig(const Config& config)
    : Config(config) {
  phiSliceWidth = 2 * std::numbers::pi_v<float> / config.nMaxPhiSlice;
}

GraphBasedTrackSeeder::GraphBasedTrackSeeder(
    const DerivedConfig& config, std::shared_ptr<GbtsGeometry> geometry,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(config),
      m_geometry(std::move(geometry)),
      m_logger(std::move(logger)) {
  if (m_cfg.useClusterWidthCuts && m_cfg.tauLookupTable.empty()) {
    throw std::invalid_argument(
        "GraphBasedTrackSeeder: the cluster width cuts need a tau lookup "
        "table");
  }
}

GbtsNodeStorage GraphBasedTrackSeeder::makeNodeStorage() const {
  GbtsNodeStorage::Config config;
  config.useClusterWidthCuts = m_cfg.useClusterWidthCuts;
  config.maxEndcapClusterWidth = m_cfg.maxEndcapClusterWidth;
  config.moduleHalfLengthY = m_cfg.moduleHalfLengthY;
  config.moduleEdgeTolerance = m_cfg.moduleEdgeTolerance;
  config.phiSliceWidth = m_cfg.phiSliceWidth;
  config.phiIndexMargin = m_cfg.phiIndexMargin;
  config.tauLutBinWidth = m_cfg.tauLutBinWidth;

  return GbtsNodeStorage(config, m_geometry, m_cfg.tauLookupTable);
}

template <typename graph_builder_t>
void GraphBasedTrackSeeder::createSeeds(const SpacePointContainer& spacePoints,
                                        const GbtsRoiDescriptor& roi,
                                        const graph_builder_t& graphBuilder,
                                        const GbtsTrackingFilter& filter,
                                        const Options& options,
                                        SeedContainer& outputSeeds) const {
  GbtsNodeStorage nodeStorage = makeNodeStorage();

  const auto layerColumn = spacePoints.column<GbtsLayerIndex>("gbtsLayerIndex");
  const auto clusterWidthColumn = spacePoints.column<float>("clusterWidth");
  const auto localPositionColumn = spacePoints.column<float>("localPositionY");

  nodeStorage.extend(spacePoints, layerColumn, clusterWidthColumn,
                     localPositionColumn);

  nodeStorage.finalize();

  createSeeds(nodeStorage, roi, graphBuilder, filter, options, outputSeeds);
}

template <typename graph_builder_t>
void GraphBasedTrackSeeder::createSeeds(
    GbtsNodeStorage& nodeStorage, const GbtsRoiDescriptor& roi,
    const graph_builder_t& graphBuilder, const GbtsTrackingFilter& filter,
    const Options& options, SeedContainer& outputSeeds,
    std::vector<std::optional<GbtsSeedFit>>* outputFits) const {
  ACTS_DEBUG("Loaded " << nodeStorage.numberOfNodes() << " graph nodes");

  // the edge the builder makes
  using edge_t = typename graph_builder_t::EdgeType;

  detail::GbtsGraph<edge_t> graph =
      graphBuilder.buildTheGraph(roi, nodeStorage, options.bFieldInZ);

  ACTS_DEBUG("Created graph with " << graph.edgeStorage.size() << " edges and "
                                   << graph.neighbours.size()
                                   << " edge links from " << graph.nDoublets
                                   << " doublets");

  if (graph.edgeStorage.empty() || graph.neighbours.empty()) {
    ACTS_WARNING("Missing edges or edge connections");
  }

  SeedCandidates candidates = findSeedCandidates(nodeStorage.nodeView(), graph,
                                                 filter, outputFits != nullptr);
  removeClones(candidates, nodeStorage.numberOfNodes());

  const std::size_t nBefore = outputSeeds.size();
  writeSeeds(candidates, nodeStorage.m_nodes.copiedFromIndexColumn().data(),
             outputSeeds, outputFits);

  ACTS_DEBUG("GBTS created " << outputSeeds.size() - nBefore << " seeds");
  if (outputSeeds.size() == nBefore) {
    ACTS_WARNING("No Seed Candidates");
  }
}

template <typename edge_t>
GraphBasedTrackSeeder::SeedCandidates GraphBasedTrackSeeder::findSeedCandidates(
    const detail::GbtsNodeView& nodeView,
    const detail::GbtsGraph<edge_t>& graph, const GbtsTrackingFilter& filter,
    const bool withFits) const {
  SeedCandidates result;
  std::vector<SeedCandidate>& vSeedCandidates = result.candidates;

  std::uint32_t seedCounter = 0;

  // Turn the fitted path of every head edge into a seed candidate.
  for (auto& [rs, edges] : extractBestPathPerEdge(nodeView, graph, filter)) {
    const edge_t& head = graph.edgeStorage[edges.front()];
    const float seedAbsEta = std::abs(-std::log(detail::edgeExpEta(head)));

    std::vector<SpacePointIndex> vN;
    vN.reserve(edges.size() + 1);
    vN.push_back(graph.edgeStorage[edges.back()].n1);
    for (const std::uint32_t e : edges | std::views::reverse) {
      vN.push_back(graph.edgeStorage[e].n2);
    }

    // a triplet is accepted if it makes it up to this point
    if (vN.size() < 3) {
      continue;
    }

    const auto origSeedSize = static_cast<std::uint32_t>(vN.size());

    const float origSeedQuality = -rs.j / origSeedSize;

    bool seedSplitFlag = (seedAbsEta < m_cfg.maxSeedSplitEta) &&
                         (origSeedSize >= m_cfg.minSplitSeedSize) &&
                         (origSeedSize <= m_cfg.maxSplitSeedSize);

    // split the seed by dropping spacepoints
    if (seedSplitFlag) {
      // 2. "drop-outs" and the original seed candidate
      std::array<std::array<SpacePointIndex, 3>, 3> triplets{};

      // triplet parameter estimate
      std::array<float, 3> invRads{};

      triplets[0] = {vN[0], vN[origSeedSize / 2], vN[origSeedSize - 1]};

      // all but the first one
      const std::vector<SpacePointIndex> dropOut1(vN.begin() + 1, vN.end());

      triplets[1] = {dropOut1[0], dropOut1[(origSeedSize - 1) / 2],
                     dropOut1[origSeedSize - 2]};

      std::vector<SpacePointIndex> dopOut2;

      dopOut2.reserve(origSeedSize - 1);

      for (std::uint32_t k = 0; k < origSeedSize; k++) {
        if (k == origSeedSize / 2) {
          continue;  // drop the middle SP in the original seed
        }

        dopOut2.emplace_back(vN[k]);
      }

      triplets[2] = {dopOut2[0], dopOut2[(origSeedSize - 1) / 2],
                     dopOut2[origSeedSize - 2]};

      for (std::uint32_t k = 0; k < invRads.size(); k++) {
        invRads[k] = estimateCurvature(nodeView, triplets[k]);
      }

      const std::array<float, 3> diffs = {std::abs(invRads[1] - invRads[0]),
                                          std::abs(invRads[2] - invRads[0]),
                                          std::abs(invRads[2] - invRads[1])};

      const bool confirmed = diffs[0] < m_cfg.maxInvRadDiff &&
                             diffs[1] < m_cfg.maxInvRadDiff &&
                             diffs[2] < m_cfg.maxInvRadDiff;

      if (confirmed) {
        seedSplitFlag = false;  // reset the flag
      }
    }

    SeedCandidate& candidate = vSeedCandidates.emplace_back();
    candidate.quality = origSeedQuality;
    candidate.cotTheta = rs.y[1];
    candidate.forSeedSplitting = seedSplitFlag;

    // The filter state is a parabola in a fixed frame, x' along the first
    // edge (inwards), and z linear in r, taken at the last node it was
    // updated with: the innermost one for a chain searched inwards, the
    // outermost strip for one whose strips were picked up outwards
    // (`pixelsInwardStripsOutward`). Taken back along the fit to the innermost
    // node, which leaves an inward chain's state as it is.
    if (withFits) {
      const detail::GbtsNodeProxy inner = nodeView[vN.front()];
      const double refX = static_cast<double>(inner.x()) * rs.c +
                          static_cast<double>(inner.y()) * rs.s;
      const double a = refX - rs.refX;
      const double x0 = rs.x[0] + rs.x[1] * a + 0.5 * rs.x[2] * a * a;
      const double x1 = rs.x[1] + rs.x[2] * a;
      const double z0 = rs.y[0] + rs.y[1] * (inner.r() - rs.refY);
      const double slope2 = 1. + x1 * x1;
      candidate.fit =
          GbtsSeedFit{.position = Vector3(refX * rs.c - x0 * rs.s,
                                          refX * rs.s + x0 * rs.c, z0),
                      // outwards: along -x', with the slope of the state
                      .phi = std::atan2(-rs.s - rs.c * x1, -rs.c + rs.s * x1),
                      .cotTheta = rs.y[1],
                      .curvature = rs.x[2] / (slope2 * std::sqrt(slope2)),
                      .pathScore = rs.j};
    }

    candidate.nodes = std::move(vN);
    candidate.edges = std::move(edges);

    result.ranking.emplace_back(origSeedQuality, seedCounter);

    ++seedCounter;
  }

  std::ranges::sort(result.ranking);

  return result;
}

void GraphBasedTrackSeeder::removeClones(SeedCandidates& candidates,
                                         const std::uint32_t nNodes) const {
  std::vector<SeedCandidate>& vSeedCandidates = candidates.candidates;
  const std::vector<std::pair<float, std::uint32_t>>& vArgSort =
      candidates.ranking;

  // hit to track associations, indexed by graph node index
  std::vector<std::uint32_t> h2t(nNodes + 1, 0);

  std::uint32_t trackId = 0;

  for (const auto& args : vArgSort) {
    const auto& seed = vSeedCandidates[args.second];
    ++trackId;

    // loop over space points indices
    for (const SpacePointIndex node : seed.nodes) {
      const std::uint32_t hitId = node + 1;

      const std::uint32_t tid = h2t[hitId];

      // unused hit or used by a lesser track
      if (tid == 0 || tid > trackId) {
        // overwrite
        h2t[hitId] = trackId;
      }
    }
  }

  std::uint32_t trackIdx = 0;

  for (const auto& args : vArgSort) {
    const auto& seed = vSeedCandidates[args.second].nodes;

    const auto nTotal = static_cast<std::uint32_t>(seed.size());

    std::uint32_t nOther = 0;

    trackId = trackIdx + 1;

    ++trackIdx;

    for (const SpacePointIndex node : seed) {
      const std::uint32_t hitId = node + 1;

      const std::uint32_t tid = h2t[hitId];

      // taken by a better candidate
      if (tid != trackId) {
        nOther++;
      }
    }

    if (nOther > m_cfg.hitShareThreshold * nTotal) {
      // reject
      vSeedCandidates[args.second].rejected = true;
    }
  }
}

void GraphBasedTrackSeeder::writeSeeds(
    const SeedCandidates& candidates,
    const std::span<const SpacePointIndex> spacePointIndices,
    SeedContainer& outputSeeds,
    std::vector<std::optional<GbtsSeedFit>>* outputFits) const {
  std::vector<SpacePointIndex> spacePoints;
  const auto write = [&](const float quality,
                         const std::optional<GbtsSeedFit>& fit) {
    auto newSeed = outputSeeds.createSeed();
    newSeed.assignSpacePointIndices(spacePoints);
    newSeed.quality() = quality;
    if (outputFits != nullptr) {
      outputFits->push_back(fit);
    }
  };

  for (const auto& args : candidates.ranking) {
    const SeedCandidate& seed = candidates.candidates[args.second];

    // identified as a clone of a better candidate, or cut
    if (seed.rejected) {
      continue;
    }

    const auto& vN = seed.nodes;

    if (!seed.forSeedSplitting) {
      spacePoints.clear();
      for (const SpacePointIndex node : vN) {
        spacePoints.push_back(spacePointIndices[node]);
      }
      write(seed.quality, seed.fit);
      continue;
    }

    // seed split into "drop-out" seeds: without the first and without the
    // middle space point
    const auto seedSize = static_cast<std::uint32_t>(vN.size());
    for (const std::size_t skipIdx : {std::size_t{0}, seedSize / 2ul}) {
      spacePoints.clear();
      for (std::uint32_t k = 0; k < seedSize; k++) {
        if (k != skipIdx) {
          spacePoints.push_back(spacePointIndices[vN[k]]);
        }
      }
      write(seed.quality, std::nullopt);
    }
  }
}

template <typename edge_t>
std::vector<GraphBasedTrackSeeder::FittedPath<edge_t>>
GraphBasedTrackSeeder::extractBestPathPerEdge(
    const detail::GbtsNodeView& nodeView,
    const detail::GbtsGraph<edge_t>& graph,
    const GbtsTrackingFilter& filter) const {
  const std::vector<edge_t>& edgeStorage = graph.edgeStorage;
  const auto nEdges = static_cast<std::uint32_t>(edgeStorage.size());
  const auto maxLength = static_cast<std::uint32_t>(m_cfg.maxChainLength);

  // Per edge, whether its inner (bit 0) and its outer (bit 1) node is a strip
  // node: looked up once here instead of on every pass below.
  constexpr std::uint8_t kInnerStrip = 1;
  constexpr std::uint8_t kOuterStrip = 2;
  const bool anyStrips = !nodeView.stripIndex.empty();
  std::vector<std::uint8_t> stripEnds(anyStrips ? nEdges : 0, 0);
  if (anyStrips) {
    for (std::uint32_t e = 0; e < nEdges; ++e) {
      const edge_t& edge = edgeStorage[e];
      stripEnds[e] = static_cast<std::uint8_t>(
          (nodeView.strip(edge.n1) != nullptr ? kInnerStrip : 0) |
          (nodeView.strip(edge.n2) != nullptr ? kOuterStrip : 0));
    }
  }

  // The edges that hold a neighbour, in order. Any other edge ends every chain
  // it is in, so it settles nothing below and heads no path of more than one
  // edge, which makes no seed unless the minimum length is one, see below.
  std::vector<std::uint32_t> linkedEdges;
  for (std::uint32_t e = 0; e < nEdges; ++e) {
    if (!graph.innerNeighbours(e).empty()) {
      linkedEdges.push_back(e);
    }
  }

  // The level of an edge is the length of the longest chain of edges going
  // outwards from it. The inner neighbours of an edge come after it in the
  // graph, as both builders order it, so one pass in edge order settles it.
  // Edges in chains longer than `maxLength` stay above it and are left out.
  std::vector<std::uint32_t> level(nEdges, 1);
  for (const std::uint32_t e : linkedEdges) {
    const std::uint32_t next = std::min(level[e] + 1, maxLength + 1);
    for (const std::uint32_t inner : graph.innerNeighbours(e)) {
      std::uint32_t& innerLevel = level[inner];
      innerLevel = std::max(innerLevel, next);
    }
  }

  // The most edges on a chain from an edge inwards to a root, the edge
  // included, stepping one level up at a time; zero for an edge left out.
  // The edges further in come later, so one pass backwards settles it.
  std::vector<std::uint32_t> edgesToRoot(nEdges);
  for (std::uint32_t e = 0; e < nEdges; ++e) {
    edgesToRoot[e] = level[e] <= maxLength ? 1 : 0;
  }
  for (const std::uint32_t e : linkedEdges | std::views::reverse) {
    if (level[e] > maxLength) {
      continue;
    }
    // A root has no settled edge further in exactly one level up: an edge
    // further in with a longer chain outwards elsewhere does not end the chain
    const std::uint32_t nextLevel = level[e] + 1;
    std::uint32_t length = 1;
    if (nextLevel <= maxLength) {
      for (const std::uint32_t inner : graph.innerNeighbours(e)) {
        if (level[inner] == nextLevel) {
          length = std::max(length, edgesToRoot[inner] + 1);
        }
      }
    }
    edgesToRoot[e] = length;
  }

  // The same over chains of strips alone, which the strip-only search of
  // `pixelsInwardStripsOutward` follows: it does not step into an edge ending
  // on a pixel, so such an edge further in, of the particle or not, must not
  // keep a strip chain from ending at its innermost strip edge.
  constexpr std::uint8_t kStripStrip = kInnerStrip | kOuterStrip;
  std::vector<std::uint32_t> stripEdgesToRoot(
      filter.m_cfg.pixelsInwardStripsOutward && anyStrips ? nEdges : 0, 0);
  if (!stripEdgesToRoot.empty()) {
    for (std::uint32_t e = nEdges; e-- > 0;) {
      if (level[e] > maxLength || stripEnds[e] != kStripStrip) {
        continue;
      }
      const std::uint32_t nextLevel = level[e] + 1;
      std::uint32_t length = 1;
      if (nextLevel <= maxLength) {
        for (const std::uint32_t inner : graph.innerNeighbours(e)) {
          if (level[inner] == nextLevel && stripEnds[inner] == kStripStrip) {
            length = std::max(length, stripEdgesToRoot[inner] + 1);
          }
        }
      }
      stripEdgesToRoot[e] = length;
    }
  }

  // The chain length a head's path needs: `addTriplets` accepts one edge short
  // at low eta
  const auto minLengthOf = [&](const edge_t& head) {
    std::uint32_t length = m_cfg.minSeedLevel;
    if (m_cfg.addTriplets && length > 0 &&
        std::abs(-std::log(detail::edgeExpEta(head))) <=
            m_cfg.maxAbsEtaAddTriplets) {
      --length;
    }
    return length;
  };

  // A mixed path on few pixels with too few strips: see `mixedMinStrips`
  const auto tooFewStrips = [&](const std::uint32_t nPixels,
                                const std::uint32_t nStrips) {
    return m_cfg.mixedMinStrips > 0 && nPixels > 0 && nStrips > 0 &&
           nPixels < m_cfg.mixedFewPixels && nStrips < m_cfg.mixedMinStrips;
  };

  // From every edge, search all paths inwards to a root, fitting on the way.
  // The complete path with the best score per space point is the candidate
  // of the edge.
  std::vector<FittedPath<edge_t>> candidates;
  std::uint32_t minLength = 0;
  std::vector<std::uint32_t> path;
  std::vector<std::uint32_t> bestPath;
  float bestQuality = 0;
  detail::GbtsEdgeState<edge_t> bestState;
  // with `pixelsInwardStripsOutward`, the strip chains do not enter the pixels
  bool stripsOnlyInward = false;

  // the path from the head inwards, `path`, against `mixedMinStrips`
  const auto pathTooFewStrips = [&]() {
    if (m_cfg.mixedMinStrips == 0 || !anyStrips) {
      return false;
    }
    std::uint32_t nStrips = 0;
    for (const std::uint32_t e : path) {
      nStrips += (stripEnds[e] & kOuterStrip) != 0 ? 1 : 0;
    }
    nStrips += (stripEnds[path.back()] & kInnerStrip) != 0 ? 1 : 0;
    const auto nPixels = static_cast<std::uint32_t>(path.size()) + 1 - nStrips;
    return tooFewStrips(nPixels, nStrips);
  };

  const auto search = [&](const auto& self, const std::uint32_t e,
                          const detail::GbtsEdgeState<edge_t>& state) -> void {
    path.push_back(e);
    // every step goes one level up, so the depth is bounded by `maxLength`
    assert(path.size() <= maxLength);

    // a strip-only chain ends at its innermost strip edge
    const std::vector<std::uint32_t>& toRoot =
        stripsOnlyInward ? stripEdgesToRoot : edgesToRoot;
    if (toRoot[e] == 1) {
      // a root completes the chain
      const float quality = -state.j / static_cast<float>(path.size() + 1);
      const bool maxCurvatureCut =
          m_cfg.postFitCurvatureFactor > 0 &&
          std::abs(state.x[2]) * m_cfg.postFitCurvatureFactor >
              filter.m_cfg.maxCurvature;
      if (!maxCurvatureCut && path.size() >= minLength &&
          quality < bestQuality && !pathTooFewStrips()) {
        bestQuality = quality;
        bestState = state;
        bestPath = path;
      }
    } else {
      for (const std::uint32_t inner : graph.innerNeighbours(e)) {
        if (level[inner] != level[e] + 1 ||
            path.size() + toRoot[inner] < minLength ||
            (stripsOnlyInward && (stripEnds[inner] & kInnerStrip) == 0)) {
          continue;
        }
        // A hit adds at most `addHit` to the score and at most
        // `edgesToRoot[inner]` hits are left, so no path through this edge can
        // beat the best one found if even the longest cannot. A shorter one
        // scores no better per point, as the best scores under `addHit` each.
        if (!bestPath.empty()) {
          const float maxScore =
              state.j +
              static_cast<float>(edgesToRoot[inner]) * filter.m_cfg.addHit;
          const float bestScore =
              -bestQuality *
              static_cast<float>(path.size() + edgesToRoot[inner] + 1);
          if (maxScore <= bestScore) {
            continue;
          }
        }
        detail::GbtsEdgeState<edge_t> innerState = state;
        if (filter.update(nodeView, edgeStorage[inner], innerState)) {
          self(self, inner, innerState);
        }
      }
    }

    path.pop_back();
  };

  // a head's state, from the head edge alone
  const auto startState = [&](const edge_t& head,
                              detail::GbtsEdgeState<edge_t>& state) {
    state.initialize(head, nodeView, filter.m_cfg.initialVarianceX,
                     filter.m_cfg.initialVarianceY,
                     filter.m_cfg.calibrateStrips,
                     filter.m_cfg.maxStripLengthFraction);
    return filter.update(nodeView, head, state);
  };

  if (filter.m_cfg.pixelsInwardStripsOutward && anyStrips) {
    // The links from an edge to the edges further out that end on a strip,
    // the ones the strip pickup follows: the same links as the inward search,
    // the other way, with the inner neighbours inverted.
    std::vector<std::uint32_t> outOffset(nEdges + 1, 0);
    const auto isLink = [&](const std::uint32_t e, const std::uint32_t inner) {
      return level[inner] <= maxLength && level[inner] == level[e] + 1;
    };
    for (std::uint32_t e = 0; e < nEdges; ++e) {
      if ((stripEnds[e] & kOuterStrip) == 0) {
        continue;
      }
      for (const std::uint32_t inner : graph.innerNeighbours(e)) {
        if (isLink(e, inner)) {
          ++outOffset[inner + 1];
        }
      }
    }
    for (std::uint32_t e = 0; e < nEdges; ++e) {
      outOffset[e + 1] += outOffset[e];
    }
    std::vector<std::uint32_t> outLinks(outOffset[nEdges]);
    std::vector<std::uint32_t> outFill(outOffset.begin(), outOffset.end() - 1);
    for (std::uint32_t e = 0; e < nEdges; ++e) {
      if ((stripEnds[e] & kOuterStrip) == 0) {
        continue;
      }
      for (const std::uint32_t inner : graph.innerNeighbours(e)) {
        if (isLink(e, inner)) {
          outLinks[outFill[inner]++] = e;
        }
      }
    }

    // the best candidate per head edge so far
    constexpr auto kNoCandidate = std::numeric_limits<std::uint32_t>::max();
    std::vector<std::uint32_t> candidateOfHead(nEdges, kNoCandidate);
    std::vector<float> qualityOfHead(nEdges,
                                     std::numeric_limits<float>::max());
    std::vector<FittedPath<edge_t>> headCandidates;
    // `edges` from the head inwards
    const auto store = [&](const std::uint32_t head,
                           const detail::GbtsEdgeState<edge_t>& state,
                           const float quality,
                           const std::vector<std::uint32_t>& edges) {
      qualityOfHead[head] = quality;
      if (candidateOfHead[head] == kNoCandidate) {
        candidateOfHead[head] =
            static_cast<std::uint32_t>(headCandidates.size());
        headCandidates.emplace_back(state, edges);
      } else {
        headCandidates[candidateOfHead[head]] =
            FittedPath<edge_t>{state, edges};
      }
    };

    // The strips picked up outwards from a pixel path, `pixelPath` from its
    // head inwards: each strip edge reached ends a candidate path of its own.
    std::vector<std::uint32_t> pixelPath;
    std::vector<std::uint32_t> stripPath;
    std::vector<std::uint32_t> fullPath;
    const auto extend =
        [&](const auto& self, const std::uint32_t e,
            const detail::GbtsEdgeState<edge_t>& state) -> void {
      stripPath.push_back(e);
      const auto total =
          static_cast<std::uint32_t>(pixelPath.size() + stripPath.size());
      const auto nPixels = static_cast<std::uint32_t>(pixelPath.size()) + 1;
      const auto nStrips = static_cast<std::uint32_t>(stripPath.size());
      const bool tooFew = tooFewStrips(nPixels, nStrips);
      const bool maxCurvatureCut =
          m_cfg.postFitCurvatureFactor > 0 &&
          std::abs(state.x[2]) * m_cfg.postFitCurvatureFactor >
              filter.m_cfg.maxCurvature;
      const float quality = -state.j / static_cast<float>(total + 1);
      if (!tooFew && !maxCurvatureCut && total >= minLengthOf(edgeStorage[e]) &&
          quality < qualityOfHead[e]) {
        fullPath.assign(stripPath.rbegin(), stripPath.rend());
        fullPath.insert(fullPath.end(), pixelPath.begin(), pixelPath.end());
        store(e, state, quality, fullPath);
      }
      // at most `level[e] - 1` strips are still to come, one per edge
      if (!tooFew || nStrips + level[e] - 1 >= m_cfg.mixedMinStrips) {
        for (std::uint32_t k = outOffset[e]; k < outOffset[e + 1]; ++k) {
          const std::uint32_t outer = outLinks[k];
          detail::GbtsEdgeState<edge_t> outerState = state;
          if (filter.updateOutward(nodeView, edgeStorage[outer], outerState)) {
            self(self, outer, outerState);
          }
        }
      }
      stripPath.pop_back();
    };

    // Each pixel head's best path, kept as a seed if long enough, then
    // extended with strips from the head's outer node.
    for (std::uint32_t head = 0; head < nEdges; ++head) {
      const edge_t& headEdge = edgeStorage[head];
      if ((stripEnds[head] & kOuterStrip) != 0 || edgesToRoot[head] == 0) {
        continue;
      }
      // A head with no strip to extend to only needs a path long enough to be
      // a seed of its own; one that can extend needs any path.
      const bool canExtend = outOffset[head] != outOffset[head + 1];
      minLength = canExtend ? 1 : minLengthOf(headEdge);
      if (edgesToRoot[head] < minLength) {
        continue;
      }
      detail::GbtsEdgeState<edge_t> state;
      if (!startState(headEdge, state)) {
        continue;
      }
      bestQuality = std::numeric_limits<float>::max();
      bestPath.clear();
      search(search, head, state);
      if (bestPath.empty()) {
        continue;
      }
      if (bestPath.size() >= minLengthOf(headEdge)) {
        store(head, bestState, bestQuality, bestPath);
      }
      // too few pixels for a path that cannot reach enough strips
      const auto nChainPixels = static_cast<std::uint32_t>(bestPath.size()) + 1;
      if (m_cfg.mixedMinStrips > 0 && nChainPixels < m_cfg.mixedFewPixels &&
          level[head] - 1 < m_cfg.mixedMinStrips) {
        continue;
      }
      pixelPath = bestPath;
      for (std::uint32_t k = outOffset[head]; k < outOffset[head + 1]; ++k) {
        const std::uint32_t outer = outLinks[k];
        detail::GbtsEdgeState<edge_t> outerState = bestState;
        if (filter.updateOutward(nodeView, edgeStorage[outer], outerState)) {
          extend(extend, outer, outerState);
        }
      }
    }

    // Chains of strips alone, inwards as usual.
    stripsOnlyInward = true;
    for (std::uint32_t head = 0; head < nEdges; ++head) {
      const edge_t& headEdge = edgeStorage[head];
      if ((stripEnds[head] & (kInnerStrip | kOuterStrip)) !=
          (kInnerStrip | kOuterStrip)) {
        continue;
      }
      minLength = minLengthOf(headEdge);
      if (stripEdgesToRoot[head] == 0 || stripEdgesToRoot[head] < minLength) {
        continue;
      }
      detail::GbtsEdgeState<edge_t> state;
      if (!startState(headEdge, state)) {
        continue;
      }
      bestQuality = std::numeric_limits<float>::max();
      bestPath.clear();
      search(search, head, state);
      if (!bestPath.empty() && bestQuality < qualityOfHead[head]) {
        store(head, bestState, bestQuality, bestPath);
      }
    }
    stripsOnlyInward = false;

    // in the order of the heads, as the inward search gives them
    candidates.reserve(headCandidates.size());
    for (std::uint32_t head = 0; head < nEdges; ++head) {
      if (candidateOfHead[head] != kNoCandidate) {
        candidates.push_back(std::move(headCandidates[candidateOfHead[head]]));
      }
    }
    return candidates;
  }

  // An edge without an inner neighbour is a path of one edge, a seed only
  // where the minimum length comes down to one: at a `minSeedLevel` of two
  // with `addTriplets`, as the triplet graph can be run.
  const bool singleEdgeSeeds =
      m_cfg.minSeedLevel <= 1 || (m_cfg.addTriplets && m_cfg.minSeedLevel <= 2);
  std::vector<std::uint32_t> allEdges;
  if (singleEdgeSeeds) {
    allEdges.resize(nEdges);
    std::iota(allEdges.begin(), allEdges.end(), 0u);
  }

  for (const std::uint32_t head : singleEdgeSeeds ? allEdges : linkedEdges) {
    const edge_t& headEdge = edgeStorage[head];

    minLength = minLengthOf(headEdge);

    if (edgesToRoot[head] == 0 || edgesToRoot[head] < minLength) {
      continue;
    }

    detail::GbtsEdgeState<edge_t> state;
    if (!startState(headEdge, state)) {
      continue;
    }

    bestQuality = std::numeric_limits<float>::max();
    bestPath.clear();
    search(search, head, state);
    if (bestPath.empty()) {
      continue;
    }

    candidates.emplace_back(bestState, bestPath);
  }

  return candidates;
}

// The seeder is written once and built for each graph builder and each edge
// that exists.
template void GraphBasedTrackSeeder::createSeeds<GbtsGraphBuilder>(
    const SpacePointContainer&, const GbtsRoiDescriptor&,
    const GbtsGraphBuilder&, const GbtsTrackingFilter&, const Options&,
    SeedContainer&) const;

template void GraphBasedTrackSeeder::createSeeds<GbtsGraphBuilder>(
    GbtsNodeStorage&, const GbtsRoiDescriptor&, const GbtsGraphBuilder&,
    const GbtsTrackingFilter&, const Options&, SeedContainer&,
    std::vector<std::optional<GbtsSeedFit>>*) const;

template GraphBasedTrackSeeder::SeedCandidates
GraphBasedTrackSeeder::findSeedCandidates<detail::GbtsEdge>(
    const detail::GbtsNodeView&, const detail::GbtsGraph<detail::GbtsEdge>&,
    const GbtsTrackingFilter&, bool) const;

template GraphBasedTrackSeeder::SeedCandidates
GraphBasedTrackSeeder::findSeedCandidates<detail::TripletGbtsEdge>(
    const detail::GbtsNodeView&,
    const detail::GbtsGraph<detail::TripletGbtsEdge>&,
    const GbtsTrackingFilter&, bool) const;

float GraphBasedTrackSeeder::estimateCurvature(
    const detail::GbtsNodeView& nodeView,
    const std::array<SpacePointIndex, 3>& nodes) const {
  // conformal mapping with the center at the last spacepoint

  std::array<float, 2> u{};
  std::array<float, 2> v{};

  const detail::GbtsNodeProxy n0 = nodeView[nodes[2]];

  const float x0 = n0.x();
  const float y0 = n0.y();

  const float r0 = n0.r();

  const float cosA = x0 / r0;

  const float sinA = y0 / r0;

  for (std::uint32_t k = 0; k < 2; k++) {
    const detail::GbtsNodeProxy nk = nodeView[nodes[k]];

    const float dx = nk.x() - x0;

    const float dy = nk.y() - y0;

    const float r2Inv = 1.0 / (dx * dx + dy * dy);

    const float xn = dx * cosA + dy * sinA;

    const float yn = -dx * sinA + dy * cosA;

    u[k] = xn * r2Inv;
    v[k] = yn * r2Inv;
  }

  const float du = u[0] - u[1];

  if (du == 0.0) {
    return 0.0;
  }

  const float A = (v[0] - v[1]) / du;

  const float B = v[1] - A * u[1];

  // curavture in units of 1/mm
  return B / std::sqrt(1 + A * A);
}

}  // namespace Acts::Experimental
