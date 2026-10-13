// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsGraphConfig.hpp"
#include "Acts/Seeding/GbtsNodeStorage.hpp"
#include "Acts/Seeding/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding/detail/DisplacedGraphTypes.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <memory>
#include <optional>

namespace Acts::Experimental {

/// The doublet graph of the GBTS workflow, built without assuming the track
/// came from the beamline.
///
/// Does the same job as @c PromptGbtsGraphBuilder and takes the same
/// configuration, but a doublet here cannot be read for a curvature, an
/// azimuth at the perigee or a @c z0: all three of those come from treating
/// the beamline as a third point on the circle, which is exactly the
/// assumption large radius tracking drops. So the doublet stage keeps to the
/// non-bending plane, but for a floor on the curvature that @c d0Max leaves
/// it, and everything else waits for a triplet fit, which is also where a
/// strip end is resolved along its strip, the fitted tangent being the first
/// direction good enough to do it with.
///
/// A doublet therefore carries the triplets it is the inner doublet of rather
/// than parameters of its own, and two doublets are matched through those.
/// Dropping the beamline also makes for many more doublets than triplets, so
/// a doublet is only held until the bins that link to its own are built, and
/// only becomes an edge of the graph once a triplet passes through it. The
/// levels and chain heads it produces are the prompt graph's, so the rest of
/// the workflow is unchanged.
class DisplacedGbtsGraphBuilder {
 public:
  /// The edge this graph builds, so that the seeder and the filter can follow
  /// a graph to its edge type rather than being told both.
  using EdgeType = detail::DisplacedGbtsEdge;

  /// Configuration, shared with the prompt graph so that a caller can hand
  /// the same object to either.
  using Config = GbtsGraphConfig;

  /// Takes the same objects as the prompt graph, differing only in the
  /// configuration it is given.
  /// @param config Configuration for the graph
  /// @param geometry GBTS geometry
  /// @param logger Logging instance
  DisplacedGbtsGraphBuilder(
      const Config& config, std::shared_ptr<const GbtsGeometry> geometry,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "DisplacedGbtsGraphBuilder", Acts::Logging::Level::INFO));

  /// Access the configuration.
  /// @return The configuration
  const Config& config() const { return m_cfg; }

  /// Build the displaced doublet graph from nodes.
  /// @param roi Region of interest descriptor
  /// @param nodeStorage Data storage containing nodes
  /// @param bFieldInZ Magnetic field in z, in GeV/(e*mm)
  /// @return The graph, with its edges and their link counts. The edges are
  ///         the doublets some triplet passed through, in no particular
  ///         order.
  detail::GbtsGraph<EdgeType> buildTheGraph(const GbtsRoiDescriptor& roi,
                                            GbtsNodeStorage& nodeStorage,
                                            float bFieldInZ) const;

 private:
  /// Fit the circle through the three nodes of a triplet.
  ///
  /// The transverse plane is inverted about the middle node, which maps every
  /// circle through that node onto a straight line, so the fit is the line
  /// through the two remaining nodes and needs no iteration.
  ///
  /// @param points The three node positions, inside out
  /// @return The circle, or nothing if the nodes are degenerate
  std::optional<detail::TripletCircle> fitTripletCircle(
      const std::array<std::array<float, 3>, 3>& points) const;

  Config m_cfg;

  std::shared_ptr<const GbtsGeometry> m_geometry;

  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace Acts::Experimental
