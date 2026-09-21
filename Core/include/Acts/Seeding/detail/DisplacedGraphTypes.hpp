
#pragma once

#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace Acts::Experimental::detail{

  /// The three fit parameters GbtsEdge kept in its @c p array, in the same
  /// order, so that slot n here means what p[n] meant there. What differs is
  /// where they come from: a displaced doublet cannot fix a circle, so these
  /// belong to a triplet the edge is part of and not to the edge itself, and
  /// an edge holds one set per such triplet.
  ///
  /// They are what a later edge is matched against. Both the recorded values
  /// and the ones matched against them describe the circle at the inner node
  /// of the doublet they share, which is the only point the two triplets
  /// agree on by construction.
  struct TripletProperties{

    /// Constructor
    /// @param expEta_ exp(eta) of the triplet
    /// @param curvature_ Signed curvature of the fitted circle
    /// @param phi_ Tangent azimuth at the triplet's innermost node
    TripletProperties(float expEta_, float curvature_, float phi_)
    : expEta{expEta_}, curvature{curvature_}, phi{phi_}{};
    /// exp(eta) of the triplet, was p[0]. Along the fitted arc rather than the
    /// chord, so it is not quite the edge's own @c expEta.
    float expEta{};
    /// Signed curvature of the fitted circle, was p[1]. In the prompt graph's
    /// dphi/dr convention, i.e. half the geometric 1/R.
    float curvature{};
    /// Azimuth of the track tangent at the innermost of the triplet's three
    /// nodes, which is this edge's own inner node, was p[2]
    float phi{};
  };

  struct DisplacedGbtsEdge final {
  DisplacedGbtsEdge() = default;

  /// Constructor
  /// @param n1_ Inner node index
  /// @param n2_ Outer node index
  /// @param expEta_ exp(eta) of the edge, the first fit parameter
  /// @param n2BarrelOrder_ Pixel barrel ordinal of the outer node's layer
  DisplacedGbtsEdge(SpacePointIndex n1_, SpacePointIndex n2_, float expEta_,
           std::int32_t n2BarrelOrder_)
      : n1{n1_},
        n2{n2_},
        level{1},
        next{1},
        expEta(expEta_),
        n2BarrelOrder{n2BarrelOrder_} {

          properties.reserve(kGbtsMaxEdgeNeighbours);
          
        }

  /// Inner node of the edge
  SpacePointIndex n1{kSpacePointIndexInvalid};
  /// Outer node of the edge
  SpacePointIndex n2{kSpacePointIndexInvalid};

  std::int8_t level{-1};
  std::int8_t next{-1};

  std::uint8_t nNei{0};

  /// exp(eta) of the doublet, from the chord between its two nodes. All a
  /// displaced edge knows about itself before a triplet fit.
  float expEta{};

  /// One entry per triplet this edge is the inner edge of, at most
  /// @c kGbtsMaxEdgeNeighbours of them. Empty while the edge is the outermost
  /// pair of every chain it could start.
  std::vector<TripletProperties> properties{};

  /// Inside-out pixel barrel ordinal of the outer node's layer, -1 for the
  /// rest. It is also the only thing the innermost neighbour loop asks about
  /// the outer node's layer, so it is cached next to the fit parameters rather
  /// than chased through the node's bin.
  std::int32_t n2BarrelOrder{-1};

  std::array<std::uint32_t, kGbtsMaxEdgeNeighbours> vNei{};
};
}