// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"

#include "Acts/Seeding/DisplacedGbtsGraphBuilder.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "Acts/Seeding/PromptGbtsGraphBuilder.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace Acts::Experimental {

namespace {

/// Order the edges so that each comes before every edge it links to.
///
/// The prompt builder creates its edges in that order already. The displaced
/// one creates an edge the first time a triplet passes through it, so an outer
/// edge can come after an inner edge an earlier triplet created. Nor does its
/// cut on the chord rule out a closed loop of links: the edges of a loop, and
/// those it links to, are not placed.
template <typename edge_t>
std::vector<std::uint32_t> edgesOutsideIn(const std::vector<edge_t>& edges) {
  const auto nEdges = static_cast<std::uint32_t>(edges.size());

  std::vector<std::uint32_t> order;
  order.reserve(nEdges);

  if constexpr (std::is_same_v<edge_t, detail::GbtsEdge>) {
    for (std::uint32_t e = 0; e < nEdges; ++e) {
      order.push_back(e);
    }
  } else {
    std::vector<std::uint32_t> nOuter(nEdges, 0);
    for (const edge_t& edge : edges) {
      for (std::uint32_t k = 0; k < edge.nNei; ++k) {
        ++nOuter[edge.vNei[k]];
      }
    }
    for (std::uint32_t e = 0; e < nEdges; ++e) {
      if (nOuter[e] == 0) {
        order.push_back(e);
      }
    }
    for (std::size_t i = 0; i < order.size(); ++i) {
      const edge_t& edge = edges[order[i]];
      for (std::uint32_t k = 0; k < edge.nNei; ++k) {
        if (--nOuter[edge.vNei[k]] == 0) {
          order.push_back(edge.vNei[k]);
        }
      }
    }
  }

  return order;
}

/// |d0| and |z0| of the helix through three nodes, inside out: the circle
/// through all three in the transverse plane, and the straight line through
/// them in z against the transverse path length from the point of closest
/// approach to the beamline.
std::array<float, 2> tripletImpact(
    const detail::GbtsNodeView& nodeView,
    const std::array<SpacePointIndex, 3>& nodes) {
  std::array<std::array<float, 3>, 3> p{};
  for (std::size_t k = 0; k < 3; ++k) {
    const detail::GbtsNodeProxy node = nodeView[nodes[k]];
    p[k] = {node.x(), node.y(), node.z()};
  }

  // twice the signed area of the triangle in the transverse plane
  const float area2 =
      2.f * (p[0][0] * (p[1][1] - p[2][1]) + p[1][0] * (p[2][1] - p[0][1]) +
             p[2][0] * (p[0][1] - p[1][1]));
  const float l01 = std::hypot(p[1][0] - p[0][0], p[1][1] - p[0][1]);
  const float l12 = std::hypot(p[2][0] - p[1][0], p[2][1] - p[1][1]);

  float d0 = 0.f;
  std::array<float, 3> s{};
  if (std::abs(area2) <= 1e-6f * l01 * l12) {
    // a straight line, through the outer two ends
    const float dx = p[2][0] - p[0][0];
    const float dy = p[2][1] - p[0][1];
    const float length = std::hypot(dx, dy);
    if (length == 0.f) {
      return {std::numeric_limits<float>::max(),
              std::numeric_limits<float>::max()};
    }
    const float ux = dx / length;
    const float uy = dy / length;
    d0 = std::abs(p[0][0] * uy - p[0][1] * ux);
    for (std::size_t k = 0; k < 3; ++k) {
      s[k] = p[k][0] * ux + p[k][1] * uy;
    }
  } else {
    std::array<float, 3> r2{};
    for (std::size_t k = 0; k < 3; ++k) {
      r2[k] = p[k][0] * p[k][0] + p[k][1] * p[k][1];
    }
    const float cx =
        (r2[0] * (p[1][1] - p[2][1]) + r2[1] * (p[2][1] - p[0][1]) +
         r2[2] * (p[0][1] - p[1][1])) /
        area2;
    const float cy =
        (r2[0] * (p[2][0] - p[1][0]) + r2[1] * (p[0][0] - p[2][0]) +
         r2[2] * (p[1][0] - p[0][0])) /
        area2;
    const float radius = std::hypot(p[0][0] - cx, p[0][1] - cy);
    d0 = std::abs(std::hypot(cx, cy) - radius);
    // the point of closest approach lies on the line from the centre through
    // the beamline
    const float phi0 = std::atan2(-cy, -cx);
    for (std::size_t k = 0; k < 3; ++k) {
      float dphi = std::atan2(p[k][1] - cy, p[k][0] - cx) - phi0;
      if (dphi > std::numbers::pi_v<float>) {
        dphi -= 2 * std::numbers::pi_v<float>;
      } else if (dphi < -std::numbers::pi_v<float>) {
        dphi += 2 * std::numbers::pi_v<float>;
      }
      s[k] = radius * std::abs(dphi);
    }
  }

  const float sMean = (s[0] + s[1] + s[2]) / 3.f;
  const float zMean = (p[0][2] + p[1][2] + p[2][2]) / 3.f;
  float sxx = 0.f;
  float sxz = 0.f;
  for (std::size_t k = 0; k < 3; ++k) {
    sxx += (s[k] - sMean) * (s[k] - sMean);
    sxz += (s[k] - sMean) * (p[k][2] - zMean);
  }
  const float tau = sxx > 0.f ? sxz / sxx : 0.f;
  return {d0, std::abs(zMean - tau * sMean)};
}

}  // namespace

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
  if (m_cfg.phiSortBuckets > GbtsNodeStorage::kMaxPhiSortBuckets) {
    throw std::invalid_argument(
        "GraphBasedTrackSeeder: phiSortBuckets exceeds the maximum");
  }

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
  config.phiSortBuckets = m_cfg.phiSortBuckets;
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

  // the edge the builder makes: a doublet with its own fit parameters for the
  // prompt graph, one carrying the triplets it belongs to for the displaced
  using edge_t = typename graph_builder_t::EdgeType;

  detail::GbtsGraph<edge_t> graph =
      graphBuilder.buildTheGraph(roi, nodeStorage, options.bFieldInZ);

  ACTS_DEBUG("Created graph with " << graph.edgeStorage.size() << " edges and "
                                   << graph.nConnections << " edge links from "
                                   << graph.nDoublets << " doublets");

  if (graph.edgeStorage.empty() || graph.nConnections == 0) {
    ACTS_WARNING("Missing edges or edge connections");
  }

  std::vector<OutputSeedProperties> vOutputSeeds;
  extractSeedsFromTheGraph(nodeStorage, graph, vOutputSeeds, filter,
                           outputFits != nullptr);

  ACTS_DEBUG("GBTS created " << vOutputSeeds.size() << " seeds");
  if (vOutputSeeds.empty()) {
    ACTS_WARNING("No Seed Candidates");
  }

  // add to output seed container
  for (const auto& seed : vOutputSeeds) {
    auto newSeed = outputSeeds.createSeed();
    newSeed.assignSpacePointIndices(seed.spacePoints);
    newSeed.quality() = seed.seedQuality;
    if (outputFits != nullptr) {
      outputFits->push_back(seed.fit);
    }
  }
}

template <typename edge_t>
void GraphBasedTrackSeeder::extractSeedsFromTheGraph(
    const GbtsNodeStorage& nodeStorage, const detail::GbtsGraph<edge_t>& graph,
    std::vector<OutputSeedProperties>& vOutputSeeds,
    const GbtsTrackingFilter& filter, const bool withFits) const {
  const detail::GbtsNodeView nodeView = nodeStorage.nodeView();

  std::vector<SeedCandidateProperties> vSeedCandidates;

  std::vector<std::pair<float, std::uint32_t>> vArgSort;

  std::uint32_t seedCounter = 0;

  // Turn the fitted path of every head edge into a seed candidate.
  for (const detail::GbtsEdgeState<edge_t>& rs :
       extractBestPathPerEdge(nodeView, graph, filter)) {
    const edge_t& head = *rs.vs.front();
    const float seedAbsEta = std::abs(-std::log(detail::edgeExpEta(head)));

    std::vector<SpacePointIndex> vN;

    for (auto sIt = rs.vs.rbegin(); sIt != rs.vs.rend(); ++sIt) {
      if (sIt == rs.vs.rbegin()) {
        vN.push_back((*sIt)->n1);
      }

      vN.push_back((*sIt)->n2);
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

    vSeedCandidates.emplace_back(origSeedQuality, false, vN, seedSplitFlag);

    if (withFits) {
      // filter state at the innermost node, x' along the first edge (inwards)
      const double slope2 = 1. + static_cast<double>(rs.x[1]) * rs.x[1];
      vSeedCandidates.back().fit = GbtsSeedFit{
          .position = Vector3(rs.refX * rs.c - rs.x[0] * rs.s,
                              rs.refX * rs.s + rs.x[0] * rs.c, rs.y[0]),
          // outwards: along -x', with the slope of the state
          .phi = std::atan2(-rs.s - rs.c * rs.x[1], -rs.c + rs.s * rs.x[1]),
          .cotTheta = rs.y[1],
          .curvature = rs.x[2] / (slope2 * std::sqrt(slope2))};
    }

    vArgSort.emplace_back(origSeedQuality, seedCounter);

    ++seedCounter;
  }

  // clone removal code goes below ...

  std::ranges::sort(vArgSort);

  if (m_cfg.maxTripletsPerMiddleNode > 0) {
    capTripletsPerMiddleNode(nodeView, nodeStorage.numberOfNodes(), vArgSort,
                             vSeedCandidates);
  }

  if (m_cfg.pairwiseCloneRemoval) {
    markClonesPairwise(nodeStorage.numberOfNodes(), vArgSort, vSeedCandidates);
  } else {
    // hit to track associations, indexed by graph node index
    std::vector<std::uint32_t> h2t(nodeStorage.numberOfNodes() + 1, 0);

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
        vSeedCandidates[args.second].isClone = true;  // reject
      }
    }
  }
  vOutputSeeds.reserve(vSeedCandidates.size());

  // drop the clones and split seeds if need be

  for (const auto& args : vArgSort) {
    const auto& seed = vSeedCandidates[args.second];

    if (seed.isClone) {
      continue;  // identified as a clone of a better candidate
    }

    const auto& vN = seed.nodes;

    if (!seed.forSeedSplitting) {
      // add seed to output

      std::vector<std::uint32_t> vSpIdx;

      vSpIdx.resize(vN.size());

      for (std::uint32_t k = 0; k < vSpIdx.size(); k++) {
        vSpIdx[k] = nodeStorage.spacePointIndex(vN[k]);
      }

      vOutputSeeds.emplace_back(seed.seedQuality, vSpIdx);
      vOutputSeeds.back().fit = seed.fit;

      continue;
    }

    // seed split into "drop-out" seeds

    const auto seedSize = static_cast<std::uint32_t>(vN.size());

    const std::array<std::size_t, 2> indices2drop = {
        0, seedSize / 2ul};  // the first and the middle

    for (const auto& skipIdx : indices2drop) {
      std::vector<std::uint32_t> newSeed;

      newSeed.reserve(seedSize - 1);

      for (std::uint32_t k = 0; k < seedSize; k++) {
        if (k == skipIdx) {
          continue;
        }

        newSeed.emplace_back(nodeStorage.spacePointIndex(vN[k]));
      }

      vOutputSeeds.emplace_back(seed.seedQuality, newSeed);
    }
  }
}

template <typename edge_t>
std::vector<detail::GbtsEdgeState<edge_t>>
GraphBasedTrackSeeder::extractBestPathPerEdge(
    const detail::GbtsNodeView& nodeView,
    const detail::GbtsGraph<edge_t>& graph,
    const GbtsTrackingFilter& filter) const {
  const std::vector<edge_t>& edgeStorage = graph.edgeStorage;
  const auto nEdges = static_cast<std::uint32_t>(edgeStorage.size());
  const auto maxLength = static_cast<std::uint32_t>(m_cfg.maxChainLength);

  // The level of an edge is the length of the longest chain of edges going
  // outwards from it. `vNei` only holds edges further in, so one pass with
  // every edge before those it links to settles it. Edges in chains longer than
  // `maxLength` stay above it and are left out, as is any edge the ordering
  // could not place.
  const std::vector<std::uint32_t> order = edgesOutsideIn(edgeStorage);
  std::vector<std::uint32_t> level(nEdges, maxLength + 1);
  for (const std::uint32_t e : order) {
    level[e] = 1;
  }
  for (const std::uint32_t e : order) {
    const edge_t& edge = edgeStorage[e];
    const std::uint32_t next = std::min(level[e] + 1, maxLength + 1);
    for (std::uint32_t k = 0; k < edge.nNei; ++k) {
      std::uint32_t& innerLevel = level[edge.vNei[k]];
      innerLevel = std::max(innerLevel, next);
    }
  }

  // The most edges on a chain from an edge inwards to a root, the edge
  // included, stepping one level up at a time; zero for an edge left out.
  // The edges further in come later, so one pass backwards settles it.
  std::vector<std::uint32_t> edgesToRoot(nEdges, 0);
  for (const std::uint32_t e : std::views::reverse(order)) {
    if (level[e] > maxLength) {
      continue;
    }
    // A root has no settled edge further in exactly one level up: an edge
    // further in with a longer chain outwards elsewhere does not end the chain
    const std::uint32_t nextLevel = level[e] + 1;
    std::uint32_t length = 1;
    if (nextLevel <= maxLength) {
      const edge_t& edge = edgeStorage[e];
      for (std::uint32_t k = 0; k < edge.nNei; ++k) {
        const std::uint32_t inner = edge.vNei[k];
        if (level[inner] == nextLevel) {
          length = std::max(length, edgesToRoot[inner] + 1);
        }
      }
    }
    edgesToRoot[e] = length;
  }

  // From every edge, search all paths inwards to a root, fitting on the way.
  // The complete path with the best score per space point is the candidate
  // of the edge.
  std::vector<detail::GbtsEdgeState<edge_t>> candidates;
  std::uint32_t minLength = 0;
  std::vector<std::uint32_t> path;
  std::vector<std::uint32_t> bestPath;
  float bestQuality = 0;
  detail::GbtsEdgeState<edge_t> bestState;

  const auto search = [&](const auto& self, const std::uint32_t e,
                          const detail::GbtsEdgeState<edge_t>& state) -> void {
    path.push_back(e);
    // every step goes one level up, so the depth is bounded by `maxLength`
    assert(path.size() <= maxLength);

    if (edgesToRoot[e] == 1) {
      // a root completes the chain
      const float quality = -state.j / static_cast<float>(path.size() + 1);
      const bool maxCurvatureCut =
          m_cfg.postFitCurvatureFactor > 0 &&
          std::abs(state.x[2]) * m_cfg.postFitCurvatureFactor >
              filter.m_cfg.maxCurvature;
      if (!maxCurvatureCut && path.size() >= minLength &&
          quality < bestQuality) {
        bestQuality = quality;
        bestState = state;
        bestPath = path;
      }
    } else {
      const edge_t& edge = edgeStorage[e];
      for (std::uint32_t k = 0; k < edge.nNei; ++k) {
        const std::uint32_t inner = edge.vNei[k];
        if (level[inner] != level[e] + 1 ||
            path.size() + edgesToRoot[inner] < minLength) {
          continue;
        }
        detail::GbtsEdgeState<edge_t> innerState = state;
        if (filter.update(nodeView, edgeStorage[inner], innerState)) {
          self(self, inner, innerState);
        }
      }
    }

    path.pop_back();
  };

  for (std::uint32_t head = 0; head < nEdges; ++head) {
    const edge_t& headEdge = edgeStorage[head];

    // `addTriplets` accepts a path one edge short at low eta
    minLength = m_cfg.minSeedLevel;
    if (m_cfg.addTriplets && minLength > 0 &&
        std::abs(-std::log(detail::edgeExpEta(headEdge))) <=
            m_cfg.maxAbsEtaAddTriplets) {
      --minLength;
    }

    if (edgesToRoot[head] == 0 || edgesToRoot[head] < minLength) {
      continue;
    }

    detail::GbtsEdgeState<edge_t> state;
    state.initialize(headEdge, nodeView, filter.m_cfg.initialVarianceX,
                     filter.m_cfg.initialVarianceY);
    if (!filter.update(nodeView, headEdge, state)) {
      continue;
    }

    bestQuality = std::numeric_limits<float>::max();
    bestPath.clear();
    search(search, head, state);
    if (bestPath.empty()) {
      continue;
    }

    bestState.vs.clear();
    for (const std::uint32_t e : bestPath) {
      bestState.vs.push_back(&edgeStorage[e]);
    }
    candidates.push_back(bestState);
  }

  return candidates;
}

// The seeder is written once and built for each graph builder that exists.
template void GraphBasedTrackSeeder::createSeeds<PromptGbtsGraphBuilder>(
    const SpacePointContainer&, const GbtsRoiDescriptor&,
    const PromptGbtsGraphBuilder&, const GbtsTrackingFilter&, const Options&,
    SeedContainer&) const;

template void GraphBasedTrackSeeder::createSeeds<PromptGbtsGraphBuilder>(
    GbtsNodeStorage&, const GbtsRoiDescriptor&, const PromptGbtsGraphBuilder&,
    const GbtsTrackingFilter&, const Options&, SeedContainer&,
    std::vector<std::optional<GbtsSeedFit>>*) const;

template void GraphBasedTrackSeeder::createSeeds<DisplacedGbtsGraphBuilder>(
    const SpacePointContainer&, const GbtsRoiDescriptor&,
    const DisplacedGbtsGraphBuilder&, const GbtsTrackingFilter&, const Options&,
    SeedContainer&) const;

template void GraphBasedTrackSeeder::createSeeds<DisplacedGbtsGraphBuilder>(
    GbtsNodeStorage&, const GbtsRoiDescriptor&,
    const DisplacedGbtsGraphBuilder&, const GbtsTrackingFilter&, const Options&,
    SeedContainer&, std::vector<std::optional<GbtsSeedFit>>*) const;

void GraphBasedTrackSeeder::capTripletsPerMiddleNode(
    const detail::GbtsNodeView& nodeView, const std::uint32_t nNodes,
    std::vector<std::pair<float, std::uint32_t>>& vArgSort,
    std::vector<SeedCandidateProperties>& vSeedCandidates) const {
  // the three space point candidates, by weight
  std::vector<std::pair<float, std::uint32_t>> triplets;
  for (const auto& [quality, index] : vArgSort) {
    const std::vector<SpacePointIndex>& nodes = vSeedCandidates[index].nodes;
    if (nodes.size() != 3) {
      continue;
    }
    const auto [d0, z0] =
        tripletImpact(nodeView, {nodes[0], nodes[1], nodes[2]});
    triplets.emplace_back(
        m_cfg.tripletD0Weight * d0 + m_cfg.tripletZ0Weight * z0, index);
  }
  std::ranges::sort(triplets);

  std::vector<std::uint32_t> perMiddleNode(nNodes, 0);
  for (const auto& [weight, index] : triplets) {
    std::uint32_t& n = perMiddleNode[vSeedCandidates[index].nodes[1]];
    if (n >= m_cfg.maxTripletsPerMiddleNode) {
      vSeedCandidates[index].isClone = true;
    } else {
      ++n;
    }
  }

  std::erase_if(vArgSort, [&vSeedCandidates](const auto& args) {
    return vSeedCandidates[args.second].isClone;
  });
}

void GraphBasedTrackSeeder::markClonesPairwise(
    const std::uint32_t nNodes,
    const std::vector<std::pair<float, std::uint32_t>>& vArgSort,
    std::vector<SeedCandidateProperties>& vSeedCandidates) const {
  // the kept candidates holding each node, as lists threaded through
  // `nextOwner`
  constexpr std::uint32_t kNone = std::numeric_limits<std::uint32_t>::max();
  std::vector<std::uint32_t> firstOwner(nNodes, kNone);
  std::vector<std::uint32_t> nextOwner;
  std::vector<std::uint32_t> owner;

  // the kept candidates sharing a node with the current one, and how many
  std::vector<std::pair<std::uint32_t, std::uint32_t>> shared;

  for (const auto& args : vArgSort) {
    SeedCandidateProperties& seed = vSeedCandidates[args.second];
    const auto nTotal = static_cast<std::uint32_t>(seed.nodes.size());

    shared.clear();
    std::uint32_t nTaken = 0;
    for (const SpacePointIndex node : seed.nodes) {
      for (std::uint32_t l = firstOwner[node]; l != kNone; l = nextOwner[l]) {
        const auto it = std::ranges::find(
            shared, owner[l], &std::pair<std::uint32_t, std::uint32_t>::first);
        if (it == shared.end()) {
          shared.emplace_back(owner[l], 1);
        } else {
          ++it->second;
        }
      }
      if (firstOwner[node] != kNone) {
        ++nTaken;
      }
    }

    std::uint32_t maxShared = 0;
    bool sharesWithLonger = false;
    for (const auto& [other, n] : shared) {
      maxShared = std::max(maxShared, n);
      sharesWithLonger =
          sharesWithLonger || vSeedCandidates[other].nodes.size() >= 4;
    }

    const bool isTriplet = nTotal == 3;
    if (maxShared > m_cfg.hitShareThreshold * nTotal ||
        (m_cfg.dropFullyTakenTriplets && isTriplet && nTaken == nTotal) ||
        (m_cfg.dropTripletsSharingLonger && isTriplet && sharesWithLonger)) {
      seed.isClone = true;
      continue;
    }

    for (const SpacePointIndex node : seed.nodes) {
      nextOwner.push_back(firstOwner[node]);
      owner.push_back(args.second);
      firstOwner[node] = static_cast<std::uint32_t>(owner.size() - 1);
    }
  }
}

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
