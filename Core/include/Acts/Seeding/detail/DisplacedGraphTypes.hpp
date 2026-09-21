
#pragma once

#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"

#include <array>
#include <cmath>
#include <cstdint>

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

    TripletProperties() = default;

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

  /// The circle three nodes put a track on, from the same conformal mapping
  /// the prompt graph validates a triplet with, but kept whole rather than
  /// reduced to a verdict.
  ///
  /// Displaced, this is the first place the track parameters exist at all: a
  /// doublet fixes a circle only by borrowing the beamline as a third point,
  /// which is exactly the assumption large radius tracking drops.
  struct TripletCircle {
    /// Slope of the conformal line. Its arctangent is the track direction at
    /// the middle node, in the frame the fit rotated into.
    float slope{};
    /// Intercept of the conformal line, which is what carries the curvature.
    float intercept{};
    /// Signed curvature, in the prompt graph's dphi/dr convention -- half the
    /// geometric 1/R -- so that the cut values tuned there carry over.
    float curvature{};
    /// Signed transverse impact parameter with respect to the beamline.
    float d0{};
    /// Azimuth of the middle node: the frame the local coordinates sit in.
    float phiMid{};
    /// cot(theta) along the fitted arc rather than along the chord.
    float tau{};
    /// exp(eta) from @ref tau.
    float expEta{};
    /// 1 / @ref expEta, the form the tau ratio is taken in.
    float invExpEta{};
    /// The three nodes in that frame, inside out, the middle one at the
    /// origin.
    std::array<std::array<float, 2>, 3> local{};

    /// Azimuth of the track tangent at one of the three nodes, pointing
    /// outwards. Two triplets sharing a doublet are matched through this, so
    /// it has to describe the circle and not the frame it was fitted in.
    /// @param k Which node, inside out
    /// @return The tangent azimuth, not wrapped
    float tangentPhi(const std::uint32_t k) const {
      // The tangent to x^2 + y^2 = 2ax + 2by at (x, y) is (b - y, x - a), and
      // 2b*intercept = 1, 2a*intercept = -slope. Scaling by 2|intercept| keeps
      // the direction and takes the straight track out of being a special
      // case.
      return phiMid + std::atan2(slope + 2.0f * intercept * local[k][0],
                                 1.0f - 2.0f * intercept * local[k][1]);
    }

    /// Track direction at one of the three nodes. Not normalised: the strip
    /// calibration only reads ratios of it.
    /// @param k Which node, inside out
    /// @return The direction in global coordinates
    std::array<float, 3> direction(const std::uint32_t k) const {
      const float phi = tangentPhi(k);
      return {std::cos(phi), std::sin(phi), tau};
    }
  };

  struct DisplacedGbtsEdge final {
  DisplacedGbtsEdge() = default;

  /// Constructor
  /// @param n1_ Inner node index
  /// @param n2_ Outer node index
  /// @param expEta_ exp(eta) of the edge, the first fit parameter
  /// @param chord_ Transverse length of the edge's own chord
  /// @param cosAlpha_ Cosine of the azimuth of that chord
  /// @param sinAlpha_ Sine of the azimuth of that chord
  /// @param n2Depth_ Barrel depth of the outer node's layer
  DisplacedGbtsEdge(SpacePointIndex n1_, SpacePointIndex n2_, float expEta_,
           float chord_, float cosAlpha_, float sinAlpha_,
           std::int32_t n2Depth_)
      : n1{n1_},
        n2{n2_},
        level{1},
        next{1},
        expEta(expEta_),
        chord(chord_),
        cosAlpha(cosAlpha_),
        sinAlpha(sinAlpha_),
        n2Depth{n2Depth_} {}

  /// Inner node of the edge
  SpacePointIndex n1{kSpacePointIndexInvalid};
  /// Outer node of the edge
  SpacePointIndex n2{kSpacePointIndexInvalid};

  std::int8_t level{-1};
  std::int8_t next{-1};

  std::uint8_t nNei{0};

  /// How many of @ref properties are filled. Zero while the edge is the
  /// outermost pair of every chain it could start.
  std::uint8_t nProperties{0};

  /// exp(eta) of the doublet, from the chord between its two nodes. All a
  /// displaced edge knows about itself before a triplet fit.
  float expEta{};

  /// Transverse length of the edge's own chord.
  ///
  /// Two edges meeting at a node bound the triplet's outer chord between them,
  /// L13 <= L12 + L23, which is what lets the turn below be cut against this
  /// pair's own reach rather than the whole detector's.
  float chord{};

  /// The edge's own chord direction in the transverse plane, (dx, dy) / chord.
  ///
  /// The turn from one edge to another is asin(curvature * L13), so it falls
  /// with pT and a wide turn is not worth fitting. Holding the direction as a
  /// unit vector rather than as an angle gives that turn's sine and cosine
  /// from a dot and a cross of two cached pairs, so the innermost loop reaches
  /// the cut without trigonometry, a wrap, or a square root.
  float cosAlpha{};
  /// @copydoc cosAlpha
  float sinAlpha{};

  /// How deep the outer node's layer sits in the barrel, -1 for an endcap.
  /// It is also all the innermost neighbour loop asks about the outer node's
  /// layer, so it is cached beside the rest rather than chased through the
  /// node's bin.
  std::int32_t n2Depth{-1};

  std::array<std::uint32_t, kGbtsMaxEdgeNeighbours> vNei{};

  /// One entry per triplet this edge is the inner edge of, at most
  /// @c kGbtsMaxEdgeNeighbours of them.
  ///
  /// Held inline and counted, exactly as @ref vNei is, and for the same
  /// reason: the cap is a small constant, so a vector here would be a heap
  /// allocation for every edge in the graph and a pointer chase to read one
  /// back.
  ///
  /// Last, and on purpose. It is over half the edge and only a triplet ever
  /// reads it, so everything the innermost loop wants is packed ahead of it
  /// into a single cache line.
  std::array<TripletProperties, kGbtsMaxEdgeNeighbours> properties{};
};

/// exp(eta) of an edge, which the displaced edge keeps in a field of its own
/// because a displaced doublet has no fit parameters to put it among.
///
/// @param edge The edge
/// @return exp(eta) of the edge
inline float edgeExpEta(const DisplacedGbtsEdge& edge) {
  return edge.expEta;
}
}