// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/TripletGbtsSeedSelector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace Acts::Experimental {

namespace {

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

TripletGbtsSeedSelector::TripletGbtsSeedSelector(
    const Config& config, std::unique_ptr<const Logger> logger)
    : m_cfg(config), m_logger(std::move(logger)) {
  if (!m_cfg.tripletQualityLeverEdges.empty() &&
      (m_cfg.minTripletQualityPerLever.size() + 1 !=
           m_cfg.tripletQualityLeverEdges.size() ||
       !std::ranges::is_sorted(m_cfg.tripletQualityLeverEdges))) {
    throw std::invalid_argument(
        "TripletGbtsSeedSelector: tripletQualityLeverEdges must be sorted and "
        "have one more entry than minTripletQualityPerLever");
  }
}

void TripletGbtsSeedSelector::createSeeds(
    const GraphBasedTrackSeeder& seeder, const TripletGbtsNodes& nodes,
    const detail::GbtsGraph<detail::TripletGbtsEdge>& graph,
    const GbtsTrackingFilter& filter, SeedContainer& outputSeeds,
    std::vector<std::optional<GbtsSeedFit>>* outputFits) const {
  const detail::GbtsNodeView nodeView = nodes.nodeView();
  const std::uint32_t nNodes = nodes.size();

  if (graph.edgeStorage.empty() || graph.neighbours.empty()) {
    ACTS_WARNING("Missing edges or edge connections");
  }

  SeedCandidates candidates =
      seeder.findSeedCandidates(nodeView, graph, filter, outputFits != nullptr);

  // read off the candidates' paths now, applied after the clone removal
  const std::vector<bool> below =
      belowTripletQuality(nodeView, graph, candidates);

  if (m_cfg.maxTripletsPerMiddleNode > 0) {
    capTripletsPerMiddleNode(nodeView, nNodes, candidates);
  }

  if (m_cfg.pairwiseCloneRemoval) {
    markClonesPairwise(nNodes, seeder.config().hitShareThreshold, candidates);
  } else {
    seeder.removeClones(candidates, nNodes);
  }

  // after the clone removal, in which such a seed still takes part; a split
  // seed is written as its halves regardless
  std::uint32_t nBelowTripletQuality = 0;
  if (!below.empty()) {
    for (const auto& [quality, index] : candidates.ranking) {
      GraphBasedTrackSeeder::SeedCandidate& seed = candidates.candidates[index];
      if (!seed.rejected && below[index] && !seed.forSeedSplitting) {
        seed.rejected = true;
        ++nBelowTripletQuality;
      }
    }
  }
  ACTS_DEBUG(nBelowTripletQuality << " seeds below the triplet quality cut");

  const std::size_t nBefore = outputSeeds.size();
  seeder.writeSeeds(candidates, nodes.spacePointIndices(), outputSeeds,
                    outputFits);
  ACTS_DEBUG("Triplet GBTS created " << outputSeeds.size() - nBefore
                                     << " seeds");
}

std::vector<bool> TripletGbtsSeedSelector::belowTripletQuality(
    const detail::GbtsNodeView& nodeView,
    const detail::GbtsGraph<detail::TripletGbtsEdge>& graph,
    const SeedCandidates& candidates) const {
  const bool applyTripletQuality =
      !graph.linkQuality.empty() &&
      (m_cfg.minTripletQuality > std::numeric_limits<float>::lowest() ||
       !m_cfg.tripletQualityLeverEdges.empty() ||
       m_cfg.minMeanTripletQualityForward >
           std::numeric_limits<float>::lowest());
  if (!applyTripletQuality) {
    return {};
  }

  // the minimum quality of a seed of three space points, from its lever arm if
  // there is a table
  const auto minTripletQuality = [&](const std::vector<SpacePointIndex>& vN) {
    const std::vector<float>& edges = m_cfg.tripletQualityLeverEdges;
    if (edges.empty()) {
      return m_cfg.minTripletQuality;
    }
    const float lever = nodeView[vN.back()].r() - nodeView[vN.front()].r();
    const auto it = std::ranges::upper_bound(edges, lever);
    if (it == edges.begin() || it == edges.end()) {
      return std::numeric_limits<float>::lowest();
    }
    return m_cfg.minTripletQualityPerLever[static_cast<std::size_t>(
        it - edges.begin() - 1)];
  };

  std::vector<bool> below(candidates.candidates.size(), false);
  for (std::size_t c = 0; c < candidates.candidates.size(); ++c) {
    const GraphBasedTrackSeeder::SeedCandidate& seed = candidates.candidates[c];
    const std::vector<std::uint32_t>& edges = seed.edges;
    // the mean quality of the path's links that carry one, taken inner link
    // first
    float sum = 0.f;
    std::uint32_t n = 0;
    for (std::size_t i = edges.size() - 1; i >= 1; --i) {
      const std::uint32_t outer = edges[i - 1];
      const auto inners = graph.innerNeighbours(outer);
      const auto it = std::ranges::find(inners, edges[i]);
      const float q =
          graph.linkQuality[graph.neighbourOffsets[outer] +
                            static_cast<std::uint32_t>(it - inners.begin())];
      if (!std::isnan(q)) {
        sum += q;
        ++n;
      }
    }
    if (n > 0) {
      const float mean = sum / static_cast<float>(n);
      below[c] = seed.nodes.size() == 3
                     ? mean < minTripletQuality(seed.nodes)
                     : std::abs(std::asinh(seed.cotTheta)) >=
                               m_cfg.tripletQualityForwardMinAbsEta &&
                           mean < m_cfg.minMeanTripletQualityForward;
    }
  }
  return below;
}

void TripletGbtsSeedSelector::capTripletsPerMiddleNode(
    const detail::GbtsNodeView& nodeView, const std::uint32_t nNodes,
    SeedCandidates& candidates) const {
  std::vector<GraphBasedTrackSeeder::SeedCandidate>& vSeedCandidates =
      candidates.candidates;
  std::vector<std::pair<float, std::uint32_t>>& vArgSort = candidates.ranking;
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
      vSeedCandidates[index].rejected = true;
    } else {
      ++n;
    }
  }

  std::erase_if(vArgSort, [&vSeedCandidates](const auto& args) {
    return vSeedCandidates[args.second].rejected;
  });
}

void TripletGbtsSeedSelector::markClonesPairwise(
    const std::uint32_t nNodes, const float hitShareThreshold,
    SeedCandidates& candidates) const {
  std::vector<GraphBasedTrackSeeder::SeedCandidate>& vSeedCandidates =
      candidates.candidates;
  const std::vector<std::pair<float, std::uint32_t>>& vArgSort =
      candidates.ranking;
  // the kept candidates holding each node, as lists threaded through
  // `nextOwner`
  constexpr std::uint32_t kNone = std::numeric_limits<std::uint32_t>::max();
  std::vector<std::uint32_t> firstOwner(nNodes, kNone);
  std::vector<std::uint32_t> nextOwner;
  std::vector<std::uint32_t> owner;

  // the kept candidates sharing a node with the current one, and how many
  std::vector<std::pair<std::uint32_t, std::uint32_t>> shared;

  for (const auto& args : vArgSort) {
    GraphBasedTrackSeeder::SeedCandidate& seed = vSeedCandidates[args.second];
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
    if (maxShared > hitShareThreshold * nTotal ||
        (m_cfg.dropFullyTakenTriplets && isTriplet && nTaken == nTotal) ||
        (m_cfg.dropTripletsSharingLonger && isTriplet && sharesWithLonger)) {
      seed.rejected = true;
      continue;
    }

    for (const SpacePointIndex node : seed.nodes) {
      nextOwner.push_back(firstOwner[node]);
      owner.push_back(args.second);
      firstOwner[node] = static_cast<std::uint32_t>(owner.size() - 1);
    }
  }
}

}  // namespace Acts::Experimental
