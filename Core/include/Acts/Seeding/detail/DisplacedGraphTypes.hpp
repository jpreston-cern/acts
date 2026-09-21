
#pragma once

#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace Acts::Experimental::detail{

  /// The three fit parameters GbtsEdge kept in its @c p array, in the same
  /// order, so that slot n here means what p[n] meant there.
  struct TripletProperties{

    TripletProperties(float expEta_, float curvature_, float phi_)
    : expEta{expEta_}, curvature{curvature_}, phi{phi_}{};
    /// exp(eta) of the edge, was p[0]
    float expEta{};
    /// Signed curvature of the edge, was p[1]
    float curvature{};
    /// Azimuthal angle extrapolated to the inner node, was p[2]
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
  
  float expEta{};

  std::vector<TripletProperties> properties{};

  /// Inside-out pixel barrel ordinal of the outer node's layer, -1 for the
  /// rest. It is also the only thing the innermost neighbour loop asks about
  /// the outer node's layer, so it is cached next to the fit parameters rather
  /// than chased through the node's bin.
  std::int32_t n2BarrelOrder{-1};

  std::array<std::uint32_t, kGbtsMaxEdgeNeighbours> vNei{};
};
}