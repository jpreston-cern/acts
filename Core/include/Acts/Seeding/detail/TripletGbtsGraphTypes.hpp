// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"

#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace Acts::Experimental::detail {

/// Maximum number of inner neighbour edges linked to per triplet graph edge.
/// The prompt graph keeps every neighbour; without the beamline constraint many
/// more triplets share a doublet, so the triplet graph caps them.
static constexpr std::uint32_t kTripletGbtsMaxEdgeNeighbours = 25;

//! [triplet edge]
/// A doublet at least one triplet passed through, which is all the graph
/// keeps: the two nodes. The links are held by the graph.
struct TripletGbtsEdge final {
  TripletGbtsEdge() = default;

  /// Constructor
  /// @param n1_ Inner node index
  /// @param n2_ Outer node index
  /// @param expEta_ exp(eta) of the edge, the first fit parameter
  TripletGbtsEdge(SpacePointIndex n1_, SpacePointIndex n2_, float expEta_)
      : n1{n1_}, n2{n2_}, expEta(expEta_) {}

  /// Inner node of the edge
  SpacePointIndex n1{kSpacePointIndexInvalid};
  /// Outer node of the edge
  SpacePointIndex n2{kSpacePointIndexInvalid};

  /// exp(eta) of the doublet, from the chord between its two nodes. All a
  /// triplet graph's edge knows about itself before a triplet fit.
  float expEta{};
};
//! [triplet edge]

/// exp(eta) of an edge, which the triplet graph's edge keeps in a field of its
/// own because its doublet has no fit parameters to put it among.
///
/// @param edge The edge
/// @return exp(eta) of the edge
inline float edgeExpEta(const TripletGbtsEdge& edge) {
  return edge.expEta;
}

/// Hand the edges and their links to the graph with every edge before the
/// edges it links to, the order the seeder reads it in. Edges in a closed loop
/// of links, and those only reached through one, are left out.
/// @param graph The graph, its edges in the order they were created
/// @param links The links as they were made, outer edge first
/// @param linkQuality Quality of each link, in the order of `links`, handed
///        to the graph as `linkQuality`; empty for none
/// @return The number of edges left out
std::uint32_t orderOutsideIn(
    GbtsGraph<TripletGbtsEdge>& graph,
    const std::vector<std::pair<std::uint32_t, std::uint32_t>>& links,
    std::span<const float> linkQuality = {});

}  // namespace Acts::Experimental::detail
