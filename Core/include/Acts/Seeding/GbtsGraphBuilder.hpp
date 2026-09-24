// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsGraphConfig.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Seeding/GbtsNodeStorage.hpp"
#include "Acts/Seeding/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace Acts::Experimental {

/// Builds the doublet graph of the GBTS workflow.
///
/// Turns a finalized `GbtsNodeStorage` into a graph whose edges are doublets
/// and whose links are the doublet pairs that a triplet cut accepted.
/// `GraphBasedTrackSeeder` grows chain levels over those links with
/// `GbtsGraph::runCCA` and walks the result to produce seeds.
///
/// The phi binning is the node storage's, so the sliding windows and the
/// indexing they slide over cannot disagree.
class GbtsGraphBuilder {
 public:
  /// The edge this graph builds, so that the seeder and the filter can follow
  /// a graph to its edge type rather than being told both.
  using EdgeType = detail::GbtsEdge;

  /// Configuration, shared with the displaced graph so that a caller can hand
  /// the same object to either.
  using Config = GbtsGraphConfig;

  /// @param config Configuration for the graph
  /// @param geometry GBTS geometry
  /// @param logger Logging instance
  GbtsGraphBuilder(const Config& config,
                   std::shared_ptr<const GbtsGeometry> geometry,
                   std::unique_ptr<const Acts::Logger> logger =
                       Acts::getDefaultLogger("GbtsGraphBuilder",
                                              Acts::Logging::Level::INFO));

  /// Access the configuration.
  /// @return The configuration
  const Config& config() const { return m_cfg; }

  /// Build doublet graph from nodes.
  /// @param roi Region of interest descriptor
  /// @param nodeStorage Data storage containing nodes
  /// @param bFieldInZ Magnetic field in z, in GeV/(e*mm)
  /// @return The graph, with its edges and their link counts
  detail::GbtsGraph<EdgeType> buildTheGraph(const GbtsRoiDescriptor& roi,
                                            GbtsNodeStorage& nodeStorage,
                                            float bFieldInZ) const;

 private:
  Config m_cfg;

  std::shared_ptr<const GbtsGeometry> m_geometry;

  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace Acts::Experimental
