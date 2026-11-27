// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/StripSpacePointCalibrationDetails.hpp"
#include "Acts/Seeding/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding/CylindricalSpacePointGrid.hpp"
#include "Acts/Seeding/DoubletSeedFinder.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Seeding/TripletSeedFinder.hpp"
#include "Acts/Seeding/TripletSeeder.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <vector>

namespace Acts::Experimental {

/// Finds triplets among GBTS graph nodes with the grid triplet seeder
/// (cylindrical space point grid, doublet and triplet finders, broad seed
/// filter), for a graph builder that links the triplets rather than searching
/// doublets and triplets itself.
///
/// The search is set up once and is const afterwards: each call keeps its own
/// grid, seeder cache and filter state, so it may run on several events at
/// once.
class GbtsGridTripletSearch {
 public:
  /// Marks a node without a stereo pair in `Inputs::stripIndex`
  static constexpr std::uint32_t kNoStrip = 0xffffffffu;

  struct Config {
    /// Space point grid; its `bFieldInZ` is set per event
    CylindricalSpacePointGrid::Config grid;
    /// Doublet finders from the middle space point to the bottom and the top
    /// ones; their `experimentCuts` are set by the search itself
    DoubletSeedFinder::Config bottomDoublets;
    DoubletSeedFinder::Config topDoublets;
    TripletSeedFinder::Config triplets;
    BroadTripletSeedFilter::Config filter;

    /// Fill the stereo pairs of the strip space points into the seeder's
    /// container, for the strip calibration of the doublets and the triplet
    /// strip checks
    bool useStripDetails = true;

    /// Radial range of the middle space points per z bin of the grid; empty to
    /// take the range of the space points themselves, less the two margins
    std::vector<std::array<float, 2>> middleRRanges;
    float middleRMinMargin = 10.f;
    float middleRMaxMargin = 10.f;

    /// Scale applied to the space point variances before the triplet checks
    float varianceScale = 1.f;

    /// Azimuthal-swing doublet cut: the azimuth of a track with impact
    /// parameter up to `doubletDPhiD0Max` changes between two radii by the
    /// curvature term plus asin(d0 / rInner) - asin(d0 / rOuter), the latter
    /// held to `doubletDPhiCap`
    bool doubletDPhiCut = false;
    float doubletDPhiD0Max = 0.f;
    float doubletDPhiConst = 0.015f;
    float doubletDPhiSlope = 2e-4f;
    float doubletDPhiCap = 10.f;

    /// Keep only the seeds that are the best of at least one of their space
    /// points
    bool seedQualitySelection = false;

    /// Drop a triplet with |eta| (from its bottom to its top space point) at
    /// or above `minQualityForwardMinAbsEta` whose quality (the seed filter's
    /// weight) is below this, before it reaches the graph: forward, low
    /// quality triplets rarely make a track, and the graph then neither links
    /// them nor makes seeds of them. The lowest float switches it off.
    float minQualityForward = std::numeric_limits<float>::lowest();
    /// See `minQualityForward`.
    float minQualityForwardMinAbsEta = 1.4f;
  };

  /// The space points to search, as graph nodes
  struct Inputs {
    /// Packed (x, y, z, r) per node
    std::span<const std::array<float, 4>> positions;
    /// (varianceR, varianceZ) per node
    std::span<const std::array<float, 2>> variances;
    /// The stereo pairs, reached through `stripIndex`
    std::span<const OuterStripSpacePointCalibrationDetails> strips;
    /// Index into `strips` per node, `kNoStrip` for a node without a pair;
    /// empty when no node has one
    std::span<const std::uint32_t> stripIndex;
    /// The nodes searched, in the order they enter the grid
    std::span<const std::uint32_t> nodes;
    /// Subtracted from the node x and y positions (and the strip centres)
    /// before the search, e.g. the beam spot when the nodes were not
    /// corrected for it
    std::array<float, 2> offset{0.f, 0.f};
    /// Layer per node; only needed with `allowedPairs`
    std::span<const GbtsLayerIndex> layers;
    /// Layer pairs a doublet may connect, `nLayers` x `nLayers`, the outer
    /// layer first; empty for any pair. When set, each grid bin is split by
    /// layer and a middle layer is combined only with the bottom layers and
    /// top layers it connects to, in the grid's usual neighbour bins.
    std::span<const std::uint8_t> allowedPairs;
    std::uint32_t nLayers = 0;
  };

  /// @param config The configuration
  /// @param logger The logger
  explicit GbtsGridTripletSearch(const Config& config,
                                 std::unique_ptr<const Logger> logger =
                                     getDefaultLogger("GbtsGridTripletSearch",
                                                      Logging::Level::INFO));

  /// Find the triplets of one event.
  /// @param inputs The nodes searched
  /// @param bFieldInZ The magnetic field along z
  /// @param triplets Filled with the triplets as nodes (bottom, middle, top),
  ///        in the order the seeder makes them
  /// @param quality Filled with the quality of each triplet
  /// @param belowQuality If not null, set to the number of triplets dropped
  ///        by `Config::minQualityForward`
  void findTriplets(const Inputs& inputs, float bFieldInZ,
                    std::vector<std::array<std::uint32_t, 3>>& triplets,
                    std::vector<float>& quality,
                    std::uint32_t* belowQuality = nullptr) const;

  /// @return The configuration
  const Config& config() const { return m_cfg; }

 private:
  const Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Logger> m_logger;
  std::unique_ptr<const Logger> m_filterLogger;
  TripletSeeder m_seeder;
};

}  // namespace Acts::Experimental
