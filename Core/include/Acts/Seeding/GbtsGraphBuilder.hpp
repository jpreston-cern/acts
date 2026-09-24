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

#include <array>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace Acts::Experimental {

/// Builds the doublet graph of the GBTS workflow.
///
/// Turns a finalized `GbtsNodeStorage` into a graph whose edges are doublets
/// and whose links are the doublet pairs that a triplet cut accepted, then
/// grows chain levels over those links with a connected component analysis.
/// `GraphBasedTrackSeeder` walks the result to produce seeds.
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

  /// Access the configuration, which also carries the chain selection that
  /// seed extraction has to agree with.
  /// @return The configuration
  const Config& config() const { return m_cfg; }

  /// Build doublet graph from nodes.
  /// @param roi Region of interest descriptor
  /// @param nodeStorage Data storage containing nodes
  /// @param bFieldInZ Magnetic field in z, in GeV/(e*mm)
  /// @return The graph, with its edges and their edge and link counts
  detail::GbtsGraph<EdgeType> buildTheGraph(const GbtsRoiDescriptor& roi,
                                            GbtsNodeStorage& nodeStorage,
                                            float bFieldInZ) const;

  /// Run connected component analysis on the graph.
  /// @param graph The graph, whose edge levels are updated
  /// @return The highest chain level any edge reached
  std::uint32_t runCCA(detail::GbtsGraph<EdgeType>& graph) const;

  /// extract edges that start a chain
  /// @param graph The graph
  /// @return The edges that start chains, ordered by length of chain, empty
  ///         if no chain reached the required level
  std::vector<detail::GbtsEdge*> extractChainHeads(
      detail::GbtsGraph<EdgeType>& graph) const;

 private:
  /// Check to see if z0 of segment is within the expected z range of the
  /// beamspot
  /// @param z0BitMask Sets allowed bins of allowed z value
  /// @param z0 Estimated z0 of segments z value at beamspot
  /// @param z0HistoCoeff Scalfactor that converts z coodindate into bin index
  /// @return Whether segment is within beamspot range
  bool checkZ0BitMask(std::uint16_t z0BitMask, float z0,
                      float z0HistoCoeff) const;

  /// Check a triplet against the pT and d0 cuts.
  /// @param nodeView View of the node positions and layers
  /// @param candidateTriplet The three graph nodes
  /// @param tripletMinPt Minimum transverse momentum
  /// @param tauRatio Tau ratio of the triplet
  /// @param tauRatioCut Tau ratio cut threshold
  /// @param bFieldInZ Magnetic field in z, in GeV/(e*mm)
  /// @return Whether the triplet is accepted
  bool validateTriplet(const detail::GbtsNodeView& nodeView,
                       const std::array<SpacePointIndex, 3>& candidateTriplet,
                       float tripletMinPt, float tauRatio, float tauRatioCut,
                       float bFieldInZ) const;

  Config m_cfg;

  std::shared_ptr<const GbtsGeometry> m_geometry;

  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace Acts::Experimental
