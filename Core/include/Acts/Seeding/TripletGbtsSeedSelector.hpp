// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/Seeding/GbtsSeedFit.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"
#include "Acts/Seeding/TripletGbtsNodes.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"
#include "Acts/Seeding/detail/TripletGbtsGraphTypes.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <vector>

namespace Acts::Experimental {

/// The seeds of the triplet GBTS graph (`TripletGbtsGraphBuilder`).
///
/// The seed candidates are found as for any GBTS graph
/// (`GraphBasedTrackSeeder::findSeedCandidates`), and written by the seeder
/// too; their selection in between is the triplet graph's own. Its graph is
/// dense, its links carry the quality of the grid triplet they were made from,
/// and a three space point seed of it was already fitted and cut on, so it
/// keeps a few triplets per middle space point, calls a candidate a clone only
/// against a single better one, and cuts on the triplet quality.
class TripletGbtsSeedSelector {
 public:
  struct Config {
    /// Keep at most this many three space point candidates per middle node,
    /// those whose helix through the three passes closest to the beamline
    /// and the interaction region by `tripletD0Weight` |d0| +
    /// `tripletZ0Weight` |z0|, before the clone removal. The filter's quality
    /// barely tells a true triplet from a fake one, which this does better.
    /// Zero for no limit.
    std::uint32_t maxTripletsPerMiddleNode = 0;
    /// See `maxTripletsPerMiddleNode`.
    float tripletD0Weight = 0.01f / UnitConstants::mm;
    /// See `maxTripletsPerMiddleNode`.
    float tripletZ0Weight = 0.005f / UnitConstants::mm;
    /// Call a seed candidate a clone only when a single better candidate
    /// that was kept holds more than the seeder's `hitShareThreshold` of its
    /// space points, rather than when better candidates between them, kept or
    /// not, hold that many (the GBTS clone removal,
    /// `GraphBasedTrackSeeder::removeClones`). In the dense triplet graph the
    /// latter removes a true three space point seed whenever two unrelated
    /// better candidates take one of its space points each.
    bool pairwiseCloneRemoval = false;
    /// With `pairwiseCloneRemoval`, still drop a three space point candidate
    /// all of whose space points kept better candidates hold.
    bool dropFullyTakenTriplets = false;
    /// With `pairwiseCloneRemoval`, also drop a three space point candidate
    /// that shares any space point with a kept better candidate of four or
    /// more: most often the same track, seeded again from a piece of it.
    bool dropTripletsSharingLonger = false;
    /// Drop a seed of three space points whose triplet quality is below this,
    /// after the clone removal, in which it still takes part. Most such seeds
    /// find no track, and the low quality ones least often. The lowest float
    /// (or below) switches it off.
    float minTripletQuality = -std::numeric_limits<float>::infinity();
    /// As `minTripletQuality`, for a seed of four or more space points with
    /// |eta| (of its fit) at or above `tripletQualityForwardMinAbsEta`, on the
    /// mean quality of its triplets that carry one.
    float minMeanTripletQualityForward =
        -std::numeric_limits<float>::infinity();
    /// See `minMeanTripletQualityForward`.
    float tripletQualityForwardMinAbsEta = 1.4f;
    /// In place of `minTripletQuality`, if not empty: the minimum triplet
    /// quality of a seed of three space points per bin of its lever arm (the
    /// radial distance from its inner to its outer space point), the bin edges
    /// in `tripletQualityLeverEdges` and one value per bin in
    /// `minTripletQualityPerLever`; no cut outside the edges. A triplet with a
    /// long lever arm needs a higher quality to be worth following.
    std::vector<float> tripletQualityLeverEdges;
    /// See `tripletQualityLeverEdges`.
    std::vector<float> minTripletQualityPerLever;
  };

  /// @param config The configuration
  /// @param logger Logging instance
  explicit TripletGbtsSeedSelector(
      const Config& config,
      std::unique_ptr<const Logger> logger =
          getDefaultLogger("TripletGbtsSeedSelector", Logging::Level::INFO));

  /// Create the seeds of a triplet graph.
  /// @param seeder The seeder finding the candidates and writing the seeds
  /// @param nodes The nodes of the graph
  /// @param graph The graph, ordered outside in
  /// @param filter Tracking filter fitting the paths
  /// @param outputSeeds Container the seeds are added to
  /// @param outputFits Fit of the tracking filter to each seed written, if not
  ///        null
  void createSeeds(
      const GraphBasedTrackSeeder& seeder, const TripletGbtsNodes& nodes,
      const detail::GbtsGraph<detail::TripletGbtsEdge>& graph,
      const GbtsTrackingFilter& filter, SeedContainer& outputSeeds,
      std::vector<std::optional<GbtsSeedFit>>* outputFits = nullptr) const;

  /// @return The configuration
  const Config& config() const { return m_cfg; }

 private:
  using SeedCandidates = GraphBasedTrackSeeder::SeedCandidates;

  /// Reject the three space point candidates beyond
  /// `maxTripletsPerMiddleNode` on their middle node, and take them out of
  /// the ranking.
  void capTripletsPerMiddleNode(const detail::GbtsNodeView& nodeView,
                                std::uint32_t nNodes,
                                SeedCandidates& candidates) const;

  /// Reject the candidates that `pairwiseCloneRemoval` and the options
  /// alongside it remove.
  void markClonesPairwise(std::uint32_t nNodes, float hitShareThreshold,
                          SeedCandidates& candidates) const;

  /// Per candidate, whether it is below the triplet quality cut.
  std::vector<bool> belowTripletQuality(
      const detail::GbtsNodeView& nodeView,
      const detail::GbtsGraph<detail::TripletGbtsEdge>& graph,
      const SeedCandidates& candidates) const;

  const Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts::Experimental
