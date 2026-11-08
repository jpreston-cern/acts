// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/DisplacedGbtsGraphBuilder.hpp"

#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/SpacePointFormation/detail/StripSpacePointCalibrationImpl.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <optional>
#include <span>
#include <utility>
#include <vector>

// Every rejection in buildTheGraph is marked `[CUT]`, so `grep "\[CUT\]"`
// lists them in the order a candidate meets them.

namespace Acts::Experimental {

namespace {

/// The doublets made with one eta bin as their inner bin.
///
/// Only the bins that link to this one ever read them, so they are released
/// once the last of those is built, and the graph is left holding just the
/// doublets a triplet passed through.
struct BinDoublets {
  /// The doublets, each node's in a contiguous range that the node's
  /// `GbtsNodeEdgeInfo` gives.
  std::vector<detail::DisplacedGbtsDoublet> doublets;
  /// Every doublet's exp(eta), kept beside them so the first tau ratio cut
  /// scans 4 bytes per doublet instead of loading whole ones. Same indices.
  std::vector<float> expEta;
  /// Every doublet's `chordEtaGradient`, inner end first, for a tau ratio
  /// mode with a pull; empty otherwise. Same indices.
  std::vector<std::array<float, 2>> etaGradient;
  /// Every doublet's two ends as the stereo cut settled them, inner first,
  /// for `settledTripletFit`; empty otherwise. Same indices.
  std::vector<std::array<std::array<float, 3>, 2>> settledEnds;
  /// The triplets the doublets are the inner doublet of, each doublet's in a
  /// contiguous range that the doublet gives.
  std::vector<detail::TripletProperties> properties;
  /// Bin groups not built yet that link to this bin.
  std::uint32_t nPendingReaders{0};
  /// Whether the bin's own group has been built.
  bool built{false};

  /// Give the memory back, which clearing alone would not.
  void release() {
    std::vector<detail::DisplacedGbtsDoublet>().swap(doublets);
    std::vector<float>().swap(expEta);
    std::vector<std::array<float, 2>>().swap(etaGradient);
    std::vector<std::array<std::array<float, 3>, 2>>().swap(settledEnds);
    std::vector<detail::TripletProperties>().swap(properties);
  }
};

/// Sliding window in phi over one non-empty eta bin, holding its phi-ordered
/// nodes directly so the innermost loop does not reach through the bin.
struct SlidingWindow {
  /// doublets with the bin as their inner bin, the candidates' outer doublets
  BinDoublets* doublets{};
  /// phi-ordered nodes of the bin, including the wrap-around duplicates
  const std::pair<float, SpacePointIndex>* phiNodes{};
  /// number of entries in @c phiNodes
  std::uint32_t numPhiNodes{};
  /// sliding window position
  std::uint32_t firstIt{};
  /// window half-width
  float deltaPhi{};
  /// The layer pair this bin came from.
  const GbtsLayerConnection* connection{};
  /// Index of that pair in the geometry.
  std::uint32_t connectionIdx{};
  /// The loosest cut a triplet through this window can be held to before its
  /// strip widening, for bounding the search of the outer doublets.
  float loosestTauRatioCut{};
  /// The loosest strip widening of such a triplet.
  float loosestTauRatioCorrStrip{};
  /// The same between two doublets whose taus are both settled.
  float loosestTauRatioCorrStripSettled{};
  /// Technology of the bin's layer.
  GbtsLayerTechnology technology{};
  /// Type of the bin's layer.
  GbtsLayerType type{};
  /// bounding boxes of the bin's phi-ordered nodes, a `kStereoBlock` of them
  /// each, for skipping a block that n1's stereo pair rules out
  const std::array<float, 6>* blockBoxes{};
};

/// Phi-ordered nodes per bounding box for the doublet stereo cut at n1.
constexpr std::uint32_t kStereoBlock = 16;

/// The range of v . d over the chords d from a point to anywhere in a box,
/// and the largest that the terms of v . d can be, for the rounding slack.
///
/// @param v The vector to dot with
/// @param from The point
/// @param box The box, lowest x, y and z then highest
/// @return The least and the most v . d, and the term magnitude
std::array<float, 3> dotRange(const std::array<float, 3>& v,
                              const std::array<float, 3>& from,
                              const std::array<float, 6>& box) {
  std::array<float, 3> range{0.f, 0.f, 0.f};
  for (std::size_t i = 0; i < 3; ++i) {
    const float low = v[i] * (box[i] - from[i]);
    const float high = v[i] * (box[3 + i] - from[i]);
    range[0] += std::min(low, high);
    range[1] += std::max(low, high);
    range[2] += std::max(std::abs(low), std::abs(high));
  }
  return range;
}

/// Whether no chord from a strip node to anywhere in a box can pass the
/// doublet stereo cut at that node, `calibrateOuterStripSpacePoint` with
/// `tolerance`: whether, over the whole box, one of its two numerators
/// exceeds what the denominator allows. The ranges are exact, as the dot
/// products are linear in the chord, and each is widened by far more than
/// the rounding either computation can make, so a chord the cut would pass
/// is never in a box ruled out.
///
/// @param strip The node's stereo pair
/// @param from The node
/// @param box The box, lowest x, y and z then highest
/// @param tolerance The strip length fraction
/// @return Whether the cut fails every chord into the box
bool missesStereoPair(
    const OuterStripSpacePointCalibrationDetailsDerived& strip,
    const std::array<float, 3>& from, const std::array<float, 6>& box,
    const float tolerance) {
  constexpr float slack = 1e-4f;

  const std::array<float, 3> scale =
      dotRange(strip.innerCrossOuterHalfVector, from, box);
  const float maxScale =
      std::max(std::abs(scale[0]), std::abs(scale[1])) + slack * scale[2];

  // the least |numerator| over the box, made smaller by the slack
  const auto minAbs = [&](const std::array<float, 3>& numerator) {
    const float least = numerator[0] > 0.f   ? numerator[0]
                        : numerator[1] < 0.f ? -numerator[1]
                                             : 0.f;
    return least - slack * numerator[2];
  };

  const float limit = maxScale * tolerance;
  return minAbs(dotRange(strip.innerToOuterSeparationCrossOuterHalfVector, from,
                         box)) > limit ||
         minAbs(dotRange(strip.innerToOuterSeparationCrossInnerHalfVector, from,
                         box)) > limit;
}

/// Whether a strip node on this layer slides radially when calibrated: true
/// for endcap strips. Barrel strips slide along z and keep their transverse
/// position.
constexpr bool slidesRadially(GbtsLayerType type,
                              GbtsLayerTechnology technology) {
  return type == GbtsLayerType::Endcap &&
         technology == GbtsLayerTechnology::Strip;
}

/// The first of a node's doublets, sorted by exp(eta), that can lie within
/// a tau ratio of `ratio` of a doublet with exp(eta) `expEta`, and the
/// exp(eta) past which none can. The ends are padded a little so that the
/// exact test, which the caller still makes, decides at the boundary.
///
/// @param sortedExpEta exp(eta) of the bin's doublets
/// @param first The node's first doublet
/// @param last One past the node's last doublet
/// @param expEta exp(eta) of the doublet being matched
/// @param ratio The tau ratio threshold
/// @return The first doublet to try, and the exp(eta) to stop past
std::pair<std::uint32_t, float> tauRatioWindow(
    const std::vector<float>& sortedExpEta, const std::uint32_t first,
    const std::uint32_t last, const float expEta, const float ratio) {
  // |e / expEta - 1| <= ratio for e in expEta * [1 - ratio, 1 + ratio]
  constexpr float pad = 1e-5f;
  const float low = expEta * (1.0f - ratio) * (1.0f - pad);
  const float high = expEta * (1.0f + ratio) * (1.0f + pad);

  const auto begin = sortedExpEta.begin();
  const auto it = std::lower_bound(begin + first, begin + last, low);

  return {static_cast<std::uint32_t>(it - begin), high};
}

/// Unit axis of a strip node's outer strip, along which its space point is
/// least well measured; zero for a pixel node.
///
/// @param strip The node's stereo pair, null for a pixel node
/// @return The axis
std::array<float, 3> stripAxis(
    const OuterStripSpacePointCalibrationDetailsDerived* strip) {
  if (strip == nullptr) {
    return {0.f, 0.f, 0.f};
  }
  const std::array<float, 3>& h = strip->outerHalfVector;
  const float norm = std::sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
  if (norm == 0.f) {
    return {0.f, 0.f, 0.f};
  }
  return {h[0] / norm, h[1] / norm, h[2] / norm};
}

/// How fast the eta of the chord between two points moves per unit length
/// either end slides along an axis: (u_z - tau (c . u)) / (c sqrt(1 + tau^2))
/// for the outer end, with c the chord's transverse direction, and the
/// negative of that for the inner one.
///
/// @param inner The inner point
/// @param outer The outer point
/// @param innerAxis The axis the inner end slides along, zero for none
/// @param outerAxis The axis the outer end slides along, zero for none
/// @return The rates for the inner and the outer end
std::array<float, 2> chordEtaGradient(const std::array<float, 3>& inner,
                                      const std::array<float, 3>& outer,
                                      const std::array<float, 3>& innerAxis,
                                      const std::array<float, 3>& outerAxis) {
  const float dx = outer[0] - inner[0];
  const float dy = outer[1] - inner[1];
  const float chord = fastHypot(dx, dy);
  if (chord == 0.f) {
    return {0.f, 0.f};
  }
  const float cx = dx / chord;
  const float cy = dy / chord;
  const float tau = (outer[2] - inner[2]) / chord;
  const float scale = 1.f / (chord * fastHypot(1.f, tau));
  const auto rate = [&](const std::array<float, 3>& u) {
    return (u[2] - tau * (cx * u[0] + cy * u[1])) * scale;
  };
  return {-rate(innerAxis), rate(outerAxis)};
}

/// Squared rate at which the difference in eta of two doublets sharing a
/// node moves, summed over the three nodes sliding independently.
///
/// @param inner The inner doublet's rates, inner end first
/// @param outer The outer doublet's rates, inner end first
/// @return The squared rate
float tripletEtaGradient2(const std::array<float, 2>& inner,
                          const std::array<float, 2>& outer) {
  const float shared = outer[0] - inner[1];
  return inner[0] * inner[0] + shared * shared + outer[1] * outer[1];
}

/// exp(eta) of the chord between two points.
///
/// @param inner The inner point
/// @param outer The outer point
/// @return exp(eta), or nothing if the two share a transverse position
std::optional<float> chordExpEta(const std::array<float, 3>& inner,
                                 const std::array<float, 3>& outer) {
  const float chord = fastHypot(outer[0] - inner[0], outer[1] - inner[1]);

  if (chord == 0.0f) {
    return std::nullopt;
  }

  const float tau = (outer[2] - inner[2]) / chord;

  return fastHypot(1.0f, tau) - tau;
}

/// Transverse arc length from the perigee out to radius r, for a track of
/// impact parameter d0 and curvature in the prompt convention (1 / 2R).
/// Straight-line distance plus the same sagitta correction the triplet fit
/// uses for its tau.
///
/// @return The arc length, or nothing if the track never reaches r
std::optional<float> arcFromPerigee(const float r, const float d0,
                                    const float curvature) {
  const float chordSquare = r * r - d0 * d0;

  if (chordSquare < 0.0f) {
    return std::nullopt;
  }

  const float chord = std::sqrt(chordSquare);
  const float sagittaTerm = curvature * chord;

  return chord * (1.0f + sagittaTerm * sagittaTerm / 6.0f);
}

/// Fit the circle through the three nodes of a triplet.
///
/// The transverse plane is inverted about the middle node, which maps every
/// circle through that node onto a straight line, so the fit is the line
/// through the two remaining nodes and needs no iteration. Kept here rather
/// than as a member so that the triplet loop, which calls it for every
/// candidate, can have it inlined, which the compiler does not judge worth
/// doing for a loop that size unless told.
///
/// @param points The three node positions, inside out
/// @return The circle, or nothing if the nodes are degenerate
[[gnu::always_inline]] inline std::optional<detail::TripletCircle>
fitTripletCircle(const std::array<std::array<float, 3>, 3>& points) {
  detail::TripletCircle circle{};

  const float x0 = points[1][0];
  const float y0 = points[1][1];
  const float r0 = fastHypot(x0, y0);

  if (r0 == 0.0f) {
    return std::nullopt;
  }

  const float cosA = x0 / r0;
  const float sinA = y0 / r0;

  circle.cosMid = cosA;
  circle.sinMid = sinA;

  // conformal mapping with the center at the middle spacepoint
  std::array<float, 2> u{};
  std::array<float, 2> v{};

  for (std::uint32_t k = 0; k < 2; k++) {
    const std::uint32_t spIdx = (k == 1) ? 2 : 0;

    const float dx = points[spIdx][0] - x0;
    const float dy = points[spIdx][1] - y0;

    const float d2 = dx * dx + dy * dy;

    if (d2 == 0.0f) {
      return std::nullopt;
    }

    const float r2Inv = 1.0f / d2;

    const float xn = dx * cosA + dy * sinA;
    const float yn = -dx * sinA + dy * cosA;

    circle.local[spIdx] = {xn, yn};

    u[k] = xn * r2Inv;
    v[k] = yn * r2Inv;
  }

  const float du = u[0] - u[1];

  if (du == 0.0f) {
    return std::nullopt;
  }

  const float A = (v[0] - v[1]) / du;
  const float B = v[1] - A * u[1];

  circle.slope = A;
  circle.intercept = B;
  circle.d0 = r0 * (B * r0 - A);
  circle.curvature = B / fastHypot(1.0f, A);

  // cot(theta) per unit transverse path, which is the arc and not the chord:
  // the sagitta a displaced track turns through over three layers is small but
  // it biases tau the same way for every triplet, so taking it out keeps the
  // triplet tau comparable with the doublet one.
  const float dx13 = points[2][0] - points[0][0];
  const float dy13 = points[2][1] - points[0][1];
  const float chord13 = fastHypot(dx13, dy13);

  if (chord13 == 0.0f) {
    return std::nullopt;
  }

  const float sagittaTerm = circle.curvature * chord13;
  const float arc13 = chord13 * (1.0f + sagittaTerm * sagittaTerm / 6.0f);

  circle.tau = (points[2][2] - points[0][2]) / arc13;

  const float hypotTau = fastHypot(1.0f, circle.tau);

  circle.expEta = hypotTau - circle.tau;
  // as in the doublet case, the sum is 1 / (hypotTau - tau) and the better
  // conditioned form of it for large tau
  circle.invExpEta = hypotTau + circle.tau;

  return circle;
}

}  // namespace

DisplacedGbtsGraphBuilder::DisplacedGbtsGraphBuilder(
    const Config& config, std::shared_ptr<const GbtsGeometry> geometry,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(config),
      m_geometry(std::move(geometry)),
      m_logger(std::move(logger)) {}

detail::GbtsGraph<detail::DisplacedGbtsEdge>
DisplacedGbtsGraphBuilder::buildTheGraph(const GbtsRoiDescriptor& roi,
                                         GbtsNodeStorage& nodeStorage,
                                         const float bFieldInZ) const {
  detail::GbtsGraph<detail::DisplacedGbtsEdge> graph;
  std::vector<detail::DisplacedGbtsEdge>& edgeStorage = graph.edgeStorage;

  // z range at the outer radius, for the outer z cuts
  const float cutZMinU =
      m_cfg.minZ0 + m_cfg.maxOuterRadius * static_cast<float>(roi.dzdrMin());
  const float cutZMaxU =
      m_cfg.maxZ0 + m_cfg.maxOuterRadius * static_cast<float>(roi.dzdrMax());

  const float tripletPtMin = m_cfg.tripletPtFraction * m_cfg.minPt;

  const float ptScale = m_cfg.tuningPt / m_cfg.minPt;

  // eta dependent curvature cuts, applied to triplets (a displaced doublet has
  // no curvature)
  const float curvatureCutHighEta = m_cfg.maxCurvatureHighEta * ptScale;
  const float curvatureCutLowEta = m_cfg.maxCurvatureLowEta * ptScale;

  // for the precuts that run before the fitted tau picks one of the two
  const float curvatureCutLoosest =
      std::max(curvatureCutLowEta, curvatureCutHighEta);
  const float curvatureCutLoosestSquare =
      curvatureCutLoosest * curvatureCutLoosest;
  // the same for a doublet with an end that calibration will still slide
  const float curvatureCutSlidingSquare = curvatureCutLoosestSquare *
                                          m_cfg.curvatureScaleEndcapStrip *
                                          m_cfg.curvatureScaleEndcapStrip;

  // if doublet taus are taken with their strip ends settled
  const bool settleTau = m_cfg.settleDoubletTau && m_cfg.doubletStereoCut &&
                         nodeStorage.hasStrips();

  // if the doublet cuts read the settled ends
  const bool settledCuts = m_cfg.settledDoubletCuts && settleTau;

  // if the first triplet fit starts from the settled ends
  const bool settledFit = m_cfg.settledTripletFit && settleTau;

  // whether every doublet made so far has a settled tau, which lets a doublet
  // that has one search its partners with the tighter window
  bool allTausSettled = true;

  // per-event parts of the doublet z0 band over d0 in [0, d0Max]
  const float d0MaxSquare = m_cfg.d0Max * m_cfg.d0Max;
  const float rOut = m_cfg.maxOuterRadius;
  const float sOutMin = std::sqrt(std::max(0.f, rOut * rOut - d0MaxSquare));

  // The doublet RZ cuts on an inner end (z1, r1) and a tau. z0 = z1 - S1 * tau,
  // with S1 = sqrt(r1^2 - d0^2) the transverse path from the perigee. d0 is
  // unknown for a doublet, so z0 and z(rOut) are bands over d0 in [0, d0Max],
  // bracketed by the two ends; `s1Min` is S1 at d0Max.
  const auto passesRZ = [&](const float z1, const float r1, const float s1Min,
                            const float tau) {
    // [CUT] max |tau|
    if (std::fabs(tau) > m_cfg.maxAbsTau) {
      return false;
    }

    if (!m_cfg.doubletFilterRZ) {
      return true;
    }

    const float z0A = z1 - r1 * tau;     // d0 = 0
    const float z0B = z1 - s1Min * tau;  // d0 = d0Max

    // [CUT] doublet z0 band must overlap [minZ0, maxZ0]
    if (std::max(z0A, z0B) < m_cfg.minZ0 || std::min(z0A, z0B) > m_cfg.maxZ0) {
      return false;
    }

    const float zOutA = z1 + (rOut - r1) * tau;
    const float zOutB = z1 + (sOutMin - s1Min) * tau;

    // [CUT] doublet z(rOut) band must overlap the RoI. The two bands may pass
    // on different d0; `tripletFilterRZ` tightens this.
    return std::max(zOutA, zOutB) >= cutZMinU &&
           std::min(zOutA, zOutB) <= cutZMaxU;
  };

  // default phi half-window, from the node storage so the two cannot disagree
  const float deltaPhi0 = 0.5f * nodeStorage.m_cfg.phiSliceWidth;

  std::uint32_t nConnections = 0;

  // edges that filled their neighbour array or their triplet property store,
  // reported as possible efficiency loss
  std::uint32_t nSaturatedNeighbours = 0;
  std::uint32_t nSaturatedProperties = 0;

  // The most edges the graph may hold, scaling with the number of space
  // points by the factor maxEdgesPerSP
  const auto nMaxEdges = static_cast<std::uint32_t>(
      m_cfg.maxEdgesPerSP * static_cast<float>(nodeStorage.numberOfNodes()));

  // whether a triplet was dropped for want of room for its edges
  bool graphFull = false;

  // doublets made, held at the moment, and the most held at any moment
  std::uint32_t nDoublets = 0;
  std::size_t nHeldDoublets = 0;
  std::size_t maxHeldDoublets = 0;

  const std::span<const GbtsBinGroup> binGroups = m_geometry->binGroups();

  // Doublets are kept per inner bin, and only while a bin group still to be
  // built links to that bin. A doublet becomes a graph edge when a triplet
  // first passes through it, so the graph never holds the rest.
  std::uint32_t numBins = 0;

  for (const GbtsBinGroup& bg : binGroups) {
    numBins = std::max(numBins, bg.bin + 1);
    for (const GbtsBinLink& link : bg.links) {
      numBins = std::max(numBins, link.bin + 1);
    }
  }

  std::vector<BinDoublets> binDoublets(numBins);

  const std::span<const GbtsLayerConnection> connections =
      m_geometry->layerConnections();

  // the loosest of what any doublet kept for a bin brings as an outer pair
  struct LoosestOuter {
    float tauRatioCut{};
    float corrStrip{};
    float corrStripSettled{};
  };
  std::vector<LoosestOuter> loosestOuter(numBins);

  for (const GbtsBinGroup& bg : binGroups) {
    for (const GbtsBinLink& link : bg.links) {
      ++binDoublets[link.bin].nPendingReaders;

      const GbtsLayerConnection& pair = connections[link.connection];
      LoosestOuter& loosest = loosestOuter[bg.bin];
      loosest.tauRatioCut =
          std::max(loosest.tauRatioCut, pair.tauRatioCutOuter);
      loosest.corrStrip = std::max(loosest.corrStrip, pair.tauRatioCorrStrip);
      loosest.corrStripSettled =
          std::max(loosest.corrStripSettled, pair.tauRatioCorrStripSettled);
    }
  }

  // Release a bin's doublets if nothing is left to read them.
  const auto releaseIfDone = [&binDoublets,
                              &nHeldDoublets](const std::uint32_t binIdx) {
    BinDoublets& bin = binDoublets[binIdx];

    if (!bin.built || bin.nPendingReaders > 0) {
      return;
    }

    nHeldDoublets -= bin.doublets.size();
    bin.release();
  };

  // Called once per bin group, built or skipped: its own bin now has all of
  // its doublets, and the bins it links to have one reader fewer.
  const auto finishGroup = [&binDoublets, &releaseIfDone, &nHeldDoublets,
                            &maxHeldDoublets](const GbtsBinGroup& bg) {
    BinDoublets& bin = binDoublets[bg.bin];
    bin.built = true;

    nHeldDoublets += bin.doublets.size();
    maxHeldDoublets = std::max(maxHeldDoublets, nHeldDoublets);

    for (const GbtsBinLink& link : bg.links) {
      --binDoublets[link.bin].nPendingReaders;
      releaseIfDone(link.bin);
    }

    // a bin nothing links to, or one whose readers were all built before it
    releaseIfDone(bg.bin);
  };

  const detail::GbtsNodeView nodeView = nodeStorage.nodeView();
  const std::span<const detail::GbtsNodeParams> params =
      nodeStorage.nodeParams();
  const std::span<detail::GbtsNodeEdgeInfo> edgeInfo =
      nodeStorage.nodeEdgeInfo();

  // A node's doublets are read before its own bin is built where a bin links
  // to itself or two link to each other, and the range read then indexes a
  // bin's doublets, so nothing may be left of a graph built on these nodes
  // before.
  for (detail::GbtsNodeEdgeInfo& info : edgeInfo) {
    info.numEdges = 0;
  }

  // reused across bin groups so that the windows are allocated once
  std::vector<SlidingWindow> phiSlidingWindow;

  // reused across nodes for putting each node's doublets in exp(eta) order
  std::vector<std::uint32_t> sortOrder;
  std::vector<detail::DisplacedGbtsDoublet> sortedDoublets;
  std::vector<std::array<float, 2>> sortedEtaGradient;
  std::vector<std::array<std::array<float, 3>, 2>> sortedSettledEnds;

  // The pull modes of the tau ratio cuts, and per node the largest
  // `chordEtaGradient` of its doublets, inner end and outer, to bound the
  // window before reading any of them.
  const bool tauPull = m_cfg.tauRatioMode != GbtsTauRatioMode::Fixed;
  const bool tauFixed = m_cfg.tauRatioMode != GbtsTauRatioMode::Pull;
  const float stripRes2 =
      m_cfg.stripAlongResolution * m_cfg.stripAlongResolution;
  const float pullFloor2 = m_cfg.tauRatioPullFloor * m_cfg.tauRatioPullFloor;
  std::vector<std::array<float, 2>> maxEtaGradient;
  if (tauPull) {
    maxEtaGradient.assign(nodeStorage.numberOfNodes(), {0.f, 0.f});
  }
  // whether a candidate passes a tau ratio cut: its tau ratio against the
  // fixed tolerance, its |log| against the pull's
  const auto passesTau = [&](const float ratio, const float fixedTolerance,
                             const float logRatio, const float pullTolerance) {
    const bool fixedPass = tauFixed && ratio <= fixedTolerance;
    const bool pullPass = tauPull && std::abs(logRatio) <= pullTolerance;
    return fixedPass || pullPass;
  };
  std::vector<float> sortedExpEta;

  // if we apply calibration to the strips
  const bool calibrate = m_cfg.calibrateStrips && nodeStorage.hasStrips();

  // if doublets are cut on their strip ends' stereo pairs
  const bool stereoCut = m_cfg.doubletStereoCut && nodeStorage.hasStrips();

  // if the stereo crossings settle an endcap strip end for the turn cut
  const bool resolveSlides = calibrate && stereoCut && m_cfg.turnCutEndcapStrip;

  // the turn cut's limit where an end slides along a radial strip
  const float curvatureCutSliding =
      curvatureCutLoosest * m_cfg.curvatureScaleEndcapStrip;

  // the |tau| below which a doublet's triplets are all held to the low-eta
  // curvature limit, negative if the turn cut does not look
  const float lowEtaTauBound =
      m_cfg.curvatureFromDoubletTau
          ? m_cfg.curvatureSplitAbsTau * (1.0f - m_cfg.curvatureSplitTauMargin)
          : -1.0f;

  // Per node, the two parts of the phi window that depend on one radius
  // alone: acos(d0Max / r), the azimuth a straight track at d0Max has walked
  // by r, and sqrt(r^2 - d0Max^2), its transverse path out to r. The window
  // of a pair is then two differences.
  const bool perPairWindow = m_cfg.useEtaBinning && m_cfg.perPairPhiWindow;
  std::vector<std::array<float, 2>> nodeWindowTerms;

  if (perPairWindow) {
    nodeWindowTerms.resize(nodeStorage.numberOfNodes());

    for (std::size_t node = 0; node < nodeWindowTerms.size(); ++node) {
      const float r = params[node].r;
      nodeWindowTerms[node] = {std::acos(std::min(1.0f, m_cfg.d0Max / r)),
                               std::sqrt(std::max(0.0f, r * r - d0MaxSquare))};
    }
  }

  // Per node, the unit axis of its outer strip for the pull's eta rates,
  // derived once here rather than once per doublet and per triplet, as the
  // stereo pairs themselves are by the node storage.
  std::vector<std::array<float, 3>> nodeStripAxis;

  if (tauPull) {
    nodeStripAxis.resize(nodeStorage.numberOfNodes());

    for (std::size_t node = 0; node < nodeStripAxis.size(); ++node) {
      nodeStripAxis[node] =
          stripAxis(nodeView.strip(static_cast<SpacePointIndex>(node)));
    }
  }

  // Per outer bin, the bounding box of each `kStereoBlock` of its phi-ordered
  // nodes, at the positions the doublet stereo cut takes its chord to.
  std::vector<std::vector<std::array<float, 6>>> binBlockBoxes;

  if (stereoCut) {
    binBlockBoxes.resize(numBins);

    for (std::uint32_t bin = 0; bin < numBins; ++bin) {
      if (binDoublets[bin].nPendingReaders == 0) {
        continue;
      }

      const auto& phiNodes = nodeStorage.etaBin(bin).phiNodes;
      std::vector<std::array<float, 6>>& boxes = binBlockBoxes[bin];
      boxes.reserve((phiNodes.size() + kStereoBlock - 1) / kStereoBlock);

      for (std::size_t k = 0; k < phiNodes.size(); ++k) {
        const SpacePointIndex node = phiNodes[k].second;
        const std::array<float, 3> point{nodeView.positions[node][0],
                                         nodeView.positions[node][1],
                                         params[node].z};
        if (k % kStereoBlock == 0) {
          boxes.push_back(
              {point[0], point[1], point[2], point[0], point[1], point[2]});
        }
        std::array<float, 6>& box = boxes.back();
        for (std::size_t i = 0; i < 3; ++i) {
          box[i] = std::min(box[i], point[i]);
          box[3 + i] = std::max(box[3 + i], point[i]);
        }
      }
    }
  }

  // a phi half window held to `maxDoubletDeltaPhi`, where there is one
  const auto capDeltaPhi = [this](const float deltaPhi) {
    return m_cfg.maxDoubletDeltaPhi > 0.0f
               ? std::min(deltaPhi, m_cfg.maxDoubletDeltaPhi)
               : deltaPhi;
  };

  // The phi window of a pair of its own radii, from the two radii and each
  // one's window terms.
  const auto pairDeltaPhi = [&](const float r1, const std::array<float, 2>& t1,
                                const float r2,
                                const std::array<float, 2>& t2) {
    const bool near = std::abs(r2 - r1) < m_cfg.phiWindowSplitDeltaRadius;
    const float pairOffset =
        near ? m_cfg.phiWindowNearOffset : m_cfg.phiWindowFarOffset;
    const float pairSlope =
        near ? m_cfg.phiWindowNearSlope : m_cfg.phiWindowFarSlope;
    return capDeltaPhi(pairOffset + std::abs(t2[0] - t1[0]) +
                       std::abs(t2[1] - t1[1]) * pairSlope * ptScale);
  };

  // the window terms of a radius, as `nodeWindowTerms` holds them per node
  const auto windowTerms = [&](const float r) {
    return std::array<float, 2>{std::acos(std::min(1.0f, m_cfg.d0Max / r)),
                                std::sqrt(std::max(0.0f, r * r - d0MaxSquare))};
  };

  // nominal position of a node, before any strip calibration
  const auto nodePoint = [&nodeView](const SpacePointIndex node) {
    const std::array<float, 4>& position = nodeView.positions[node];
    return std::array<float, 3>{position[0], position[1], position[2]};
  };

  // loop over bin groups
  for (const GbtsBinGroup& bg : binGroups) {
    // B1 is the inner bin of every pair in the group
    const detail::GbtsEtaBinInfo& B1 = nodeStorage.etaBin(bg.bin);

    if (B1.empty()) {
      finishGroup(bg);
      continue;
    }

    // the doublets this group makes
    BinDoublets& inner = binDoublets[bg.bin];

    // used for the phi window width
    const float rb1 = B1.minRadius;

    const bool isPixel1 = B1.technology == GbtsLayerTechnology::Pixel;
    // barrel depth over all barrel layers, -1 for an endcap
    const std::int32_t depth1 = B1.depth;

    const bool useMatchBeforeCreate = m_cfg.matchBeforeCreate && depth1 >= 0 &&
                                      depth1 <= m_cfg.matchBeforeCreateMaxDepth;

    // prepare a sliding window for each non-empty bin2 in the group
    phiSlidingWindow.clear();

    for (const GbtsBinLink& link : bg.links) {
      const std::uint32_t b2Idx = link.bin;
      // outer bin
      const detail::GbtsEtaBinInfo& B2 = nodeStorage.etaBin(b2Idx);

      if (B2.empty()) {
        continue;
      }

      const float rb2 = B2.maxRadius;

      float deltaPhi = deltaPhi0;  // the default

      if (m_cfg.useEtaBinning) {
        const float absDr = std::fabs(rb2 - rb1);
        const float maxD0 = m_cfg.d0Max;
        auto phiWindow = [&rb1, &rb2, &maxD0, ptScale](
                             const float& phiWindowOffset,
                             const float& phiWindowSlope) {
          const float maxD0Square = maxD0 * maxD0;

          // Azimuth a straight track of impact parameter d0 walks between
          // the radii, exactly acos(d0 / r). The ratio is capped at 1 for
          // bins inside d0Max, otherwise acos gives NaN and the window breaks.
          const float frac1 = std::min(1.0f, maxD0 / rb1);
          const float frac2 = std::min(1.0f, maxD0 / rb2);
          const float displacementTerm = std::acos(frac2) - std::acos(frac1);

          const float corr1 =
              std::sqrt(std::max(0.0f, rb1 * rb1 - maxD0Square));
          const float corr2 =
              std::sqrt(std::max(0.0f, rb2 * rb2 - maxD0Square));

          const float curvatureTerm =
              (corr2 - corr1) * phiWindowSlope * ptScale;

          // displacement and bending are independent, so the worst case is
          // their sum
          const float absPhiWindow = phiWindowOffset +
                                     std::abs(displacementTerm) +
                                     std::abs(curvatureTerm);

          return absPhiWindow;
        };
        if (absDr < m_cfg.phiWindowSplitDeltaRadius) {
          deltaPhi =
              phiWindow(m_cfg.phiWindowNearOffset, m_cfg.phiWindowNearSlope);
        } else {
          deltaPhi =
              phiWindow(m_cfg.phiWindowFarOffset, m_cfg.phiWindowFarSlope);
        }
        deltaPhi = capDeltaPhi(deltaPhi);
      }

      SlidingWindow& window = phiSlidingWindow.emplace_back();
      window.doublets = &binDoublets[b2Idx];
      window.phiNodes = B2.phiNodes.data();
      window.numPhiNodes = static_cast<std::uint32_t>(B2.phiNodes.size());
      window.deltaPhi = deltaPhi;
      const GbtsLayerConnection& pair = connections[link.connection];
      const LoosestOuter& outerLoosest = loosestOuter[b2Idx];
      window.connection = &pair;
      window.connectionIdx = link.connection;
      window.loosestTauRatioCut =
          std::max(pair.tauRatioCutInner, outerLoosest.tauRatioCut);
      if (nodeStorage.hasStrips()) {
        window.loosestTauRatioCorrStripSettled = std::max(
            pair.tauRatioCorrStripSettled, outerLoosest.corrStripSettled);
        // a settled pair can be among the candidates of an unsettled search
        window.loosestTauRatioCorrStrip =
            std::max({pair.tauRatioCorrStrip, outerLoosest.corrStrip,
                      window.loosestTauRatioCorrStripSettled});
      }
      window.technology = B2.technology;
      window.type = B2.type;
      if (stereoCut) {
        window.blockBoxes = binBlockBoxes[b2Idx].data();
      }
    }

    // outer loop over n1 nodes in the inner bin
    for (SpacePointIndex n1Idx = B1.nodes.first; n1Idx < B1.nodes.second;
         ++n1Idx) {
      // initialization using the top watermark of the bin's doublets
      edgeInfo[n1Idx].firstEdge =
          static_cast<std::uint32_t>(inner.doublets.size());

      // incoming edges created for n1
      std::uint16_t numCreatedEdges = 0;

      const detail::GbtsNodeParams& n1pars = params[n1Idx];

      const float phi1 = n1pars.phi;
      const float r1 = n1pars.r;
      const float z1 = n1pars.z;

      const std::array<float, 4>& position = nodeView.positions[n1Idx];
      const float x1 = position[0];
      const float y1 = position[1];

      // this node as the triplet fit wants it, before any strip resolution
      const std::array<float, 3> point1{x1, y1, z1};

      // the stereo pair n1 was formed from, null for a pixel node
      const auto* strip1 = nodeView.strip(n1Idx);

      // n1's share of the phi window of each pair it makes
      const std::array<float, 2> windowTerms1 =
          perPairWindow ? nodeWindowTerms[n1Idx] : std::array<float, 2>{};

      // transverse path from the perigee to n1 at d0 = d0Max, for the z0 band
      const float s1Min = std::sqrt(std::max(0.f, r1 * r1 - d0MaxSquare));

      // loop over the sliding windows of the outer bins
      for (auto& slw : phiSlidingWindow) {
        const bool isPixel2 = slw.technology == GbtsLayerTechnology::Pixel;

        // Whether either end of the doublet will slide radially under
        // calibration (an endcap strip), taking the chord with it. The
        // doublet curvature cut is then loosened, and the turn cut waits for
        // the ends to settle.
        const bool doubletSlides =
            calibrate && (slidesRadially(B1.type, B1.technology) ||
                          slidesRadially(slw.type, slw.technology));

        // The turn cut is skipped when switched off, or when an end will
        // slide and the stereo crossings do not settle it. Settled here for
        // the inner two nodes; the third is checked per candidate.
        const bool slides1 = slidesRadially(B1.type, B1.technology);
        const bool slides2 = slidesRadially(slw.type, slw.technology);
        const bool innerSkipTurn =
            !m_cfg.useTurnAngleCut || (doubletSlides && !resolveSlides);

        // whether the RZ cuts wait for the stereo cut to settle the ends; two
        // pixel nodes settle to where they are, so theirs need not
        const bool settledRZ = settledCuts && !(isPixel1 && isPixel2);

        // and the transverse ones: only an endcap strip end settles across
        // the transverse plane, a barrel strip end along z
        const bool settledXY = settledCuts && doubletSlides;

        // deltaPhi is the half width: displacement and charge are both
        // symmetric about phi1
        const float deltaPhi = slw.deltaPhi;
        const float minPhi = phi1 - deltaPhi;
        const float maxPhi = phi1 + deltaPhi;

        // whether whole blocks of candidates can be ruled out by n1's stereo
        // pair, and the first candidate of the next block to try
        const bool skipBlocks = stereoCut && strip1 != nullptr;
        std::uint32_t nextBlockStart = 0;

        for (std::uint32_t n2PhiIdx = slw.firstIt; n2PhiIdx < slw.numPhiNodes;
             ++n2PhiIdx) {
          const float phi2 = slw.phiNodes[n2PhiIdx].first;

          // [CUT] phi window
          if (phi2 < minPhi) {
            // update the window position
            slw.firstIt = n2PhiIdx;
            continue;
          }
          if (phi2 > maxPhi) {
            // break and go to the next window
            break;
          }

          // [CUT] the doublet stereo cut at n1, below, for a whole block of
          // candidates at once: n1's stereo pair alone rejects most of the
          // window, and those it passes lie together in phi
          if (skipBlocks && n2PhiIdx >= nextBlockStart) {
            const std::uint32_t block = n2PhiIdx / kStereoBlock;
            nextBlockStart = (block + 1) * kStereoBlock;

            if (missesStereoPair(*strip1, point1, slw.blockBoxes[block],
                                 m_cfg.doubletStripLengthFraction)) {
              n2PhiIdx = nextBlockStart - 1;
              continue;
            }
          }

          const SpacePointIndex n2Idx = slw.phiNodes[n2PhiIdx].second;

          // [CUT] phi window of the pair's own radii, the bins' window having
          // been sized for the extremes of both bins; on the settled ends
          // below instead where an end will settle radially
          if (perPairWindow && !settledXY &&
              std::abs(phi2 - phi1) > pairDeltaPhi(r1, windowTerms1,
                                                   params[n2Idx].r,
                                                   nodeWindowTerms[n2Idx])) {
            continue;
          }

          // n2's doublets, as indices into its own bin's
          const detail::GbtsNodeEdgeInfo& n2Info = edgeInfo[n2Idx];

          const std::uint32_t n2FirstEdge = n2Info.firstEdge;
          const std::uint16_t n2NumEdges = n2Info.numEdges;
          const std::uint32_t n2LastEdge = n2FirstEdge + n2NumEdges;

          const detail::GbtsNodeParams& n2pars = params[n2Idx];

          const std::array<float, 4>& position2 = nodeView.positions[n2Idx];

          const float dx = position2[0] - x1;
          const float dy = position2[1] - y1;

          const float chord = fastHypot(dx, dy);

          // [CUT] minimum chord. Shares `minDeltaRadius` with the prompt
          // graph's radial step cut; on the settled ends below instead where
          // an end will settle radially, but never a chord of nothing.
          if (settledXY ? chord == 0.0f : chord < m_cfg.minDeltaRadius) {
            continue;
          }

          const float z2 = n2pars.z;

          // Doublet tau on nominal ends: strip calibration needs the track
          // direction, which only the triplet fit gives. `tauRatioCorrStrip`
          // covers the difference.
          const float dz = z2 - z1;
          const float tau = dz / chord;
          const float ftau = std::fabs(tau);

          // [CUT] cluster width tau window, on both nodes. Only pixel barrel
          // nodes have a finite window; the rest pass.
          const auto passesClusterWidth = [&n1pars, &n2pars](const float t) {
            return t >= n1pars.minTau && t <= n1pars.maxTau &&
                   t >= n2pars.minTau && t <= n2pars.maxTau;
          };

          // [CUT] the cluster width tau window, max |tau|, the doublet z0 and
          // z(rOut) bands, on the nominal ends unless they are to be taken on
          // the settled ones below
          if (!settledRZ &&
              (!passesClusterWidth(ftau) || !passesRZ(z1, r1, s1Min, tau))) {
            continue;
          }

          // [CUT] doublet stereo: the chord must cross both strips of the
          // stereo pair at each strip end. The triplet fit makes the same
          // check against the fitted tangent, at the tighter
          // `maxStripLengthFraction`, and it rejects nearly every candidate
          // there; taken here it spares the triplet loop those doublets.
          // Taken straight after the RZ cuts, as `DoubletSeedFinder` takes its
          // strip calibration: it removes most of what reaches it, for a few
          // dot products, and so spares the curvature cut's square root.
          // The chord the turn cut reads, which for an end sliding along a
          // radial strip is taken from where the stereo check finds the
          // chord crosses that strip rather than from the nominal end.
          float turnDx = dx;
          float turnDy = dy;
          // whether the turn cut can trust this doublet's chord
          bool turnSettled = !doubletSlides;

          // the doublet's tau, with its strip ends settled if they can be,
          // and whether they all were
          float doubletTau = tau;
          bool tauSettled = isPixel1 && isPixel2;

          // the ends that tau is taken between, nominal until settled
          std::array<float, 3> settled1{x1, y1, z1};
          std::array<float, 3> settled2{position2[0], position2[1], z2};

          if (stereoCut) {
            const std::array<float, 3> chordDirection{dx, dy, dz};
            std::array<float, 3> crossing1{};
            std::array<float, 3> crossing2{};

            if (strip1 != nullptr &&
                !Acts::detail::calibrateOuterStripSpacePoint(
                    chordDirection, *strip1, crossing1,
                    m_cfg.doubletStripLengthFraction)) {
              continue;
            }

            const auto* strip2 = nodeView.strip(n2Idx);

            if (strip2 != nullptr &&
                !Acts::detail::calibrateOuterStripSpacePoint(
                    chordDirection, *strip2, crossing2,
                    m_cfg.doubletStripLengthFraction)) {
              continue;
            }

            // a sliding end is settled only if it had a stereo pair to cross
            if (resolveSlides && doubletSlides &&
                (!slides1 || strip1 != nullptr) &&
                (!slides2 || strip2 != nullptr)) {
              const float settledX1 = slides1 ? crossing1[0] : x1;
              const float settledY1 = slides1 ? crossing1[1] : y1;
              const float settledX2 = slides2 ? crossing2[0] : position2[0];
              const float settledY2 = slides2 ? crossing2[1] : position2[1];
              turnDx = settledX2 - settledX1;
              turnDy = settledY2 - settledY1;
              turnSettled = true;
            }

            // Settle each strip end along its strip for the tau: a barrel
            // strip runs along z, an endcap strip radially. The other
            // coordinates stay nominal, as the crossing lies on the outer
            // strip of the pair rather than at the space point.
            if (settleTau && (isPixel1 || strip1 != nullptr) &&
                (isPixel2 || strip2 != nullptr)) {
              if (strip1 != nullptr) {
                if (slides1) {
                  settled1[0] = crossing1[0];
                  settled1[1] = crossing1[1];
                } else {
                  settled1[2] = crossing1[2];
                }
              }

              if (strip2 != nullptr) {
                if (slides2) {
                  settled2[0] = crossing2[0];
                  settled2[1] = crossing2[1];
                } else {
                  settled2[2] = crossing2[2];
                }
              }

              const float settledChord = fastHypot(settled2[0] - settled1[0],
                                                   settled2[1] - settled1[1]);

              if (settledChord > 0.0f) {
                doubletTau = (settled2[2] - settled1[2]) / settledChord;
                tauSettled = true;
              }
            }
          }

          // [CUT] max |tau|, the doublet z0 and z(rOut) bands, on the settled
          // tau and inner end. An endcap strip inner end has moved radially,
          // which moves its path from the perigee with it.
          if (settledRZ) {
            const bool moved1 = tauSettled && slides1 && strip1 != nullptr;
            const float r1s = moved1 ? fastHypot(settled1[0], settled1[1]) : r1;
            const float s1MinS =
                moved1 ? std::sqrt(std::max(0.f, r1s * r1s - d0MaxSquare))
                       : s1Min;

            if (!passesClusterWidth(std::fabs(doubletTau)) ||
                !passesRZ(settled1[2], r1s, s1MinS, doubletTau)) {
              continue;
            }
          }

          // [CUT] minimum chord and the pair's phi window, on the settled
          // transverse ends, which an endcap strip end has moved
          if (settledXY) {
            const float sdx = settled2[0] - settled1[0];
            const float sdy = settled2[1] - settled1[1];

            if (fastHypot(sdx, sdy) < m_cfg.minDeltaRadius) {
              continue;
            }

            if (perPairWindow) {
              const float r1s = fastHypot(settled1[0], settled1[1]);
              const float r2s = fastHypot(settled2[0], settled2[1]);
              // the azimuth between the two ends, free of any wrap
              const float dPhi = std::atan2(
                  settled1[0] * settled2[1] - settled1[1] * settled2[0],
                  settled1[0] * settled2[0] + settled1[1] * settled2[1]);

              if (std::abs(dPhi) >
                  pairDeltaPhi(r1s, windowTerms(r1s), r2s, windowTerms(r2s))) {
                continue;
              }
            }
          }

          // Two points fix no circle, but of the circles through them that
          // pass within d0Max of the beamline one bends least, and a track
          // inside the curvature limit needs that one inside it too. This is
          // the interaction point cut of `DoubletSeedFinder`, taken exactly:
          // its small d0Max / r forms do not hold where d0Max nears r1.
          //
          // A chord whose own line passes within d0Max is straight enough for
          // any limit. Otherwise, with the centre a distance t along the
          // chord's bisector, |d0| = d0Max is a quadratic in t and its larger
          // root is the widest circle. Lengths are carried times the chord, so
          // there is no division. Costs a square root, hence after the RZ and
          // stereo cuts.
          //
          // On nominal positions an endcap strip end has still to slide, which
          // turns the chord, so such a doublet is held to a limit
          // `curvatureScaleEndcapStrip` times looser instead. With
          // `settledDoubletCuts` the cut takes the settled transverse ends,
          // and the plain limit where every end settled.
          const bool curvatureSettled = settledCuts && tauSettled;
          const float cx1 = curvatureSettled ? settled1[0] : x1;
          const float cy1 = curvatureSettled ? settled1[1] : y1;
          const float cx2 = curvatureSettled ? settled2[0] : position2[0];
          const float cy2 = curvatureSettled ? settled2[1] : position2[1];
          const float cdx = cx2 - cx1;
          const float cdy = cy2 - cy1;
          const float cChord = curvatureSettled ? fastHypot(cdx, cdy) : chord;

          // chord times the distance from the beamline to the chord's line
          const float chordTimesMiss = std::abs(cx1 * cdy - cy1 * cdx);
          const float chordTimesD0Max = m_cfg.d0Max * cChord;

          if (chordTimesMiss > chordTimesD0Max) {
            const float delta = chordTimesMiss * chordTimesMiss -
                                chordTimesD0Max * chordTimesD0Max;
            const float dot = cx1 * cx2 + cy1 * cy2 - d0MaxSquare;
            // 2 * t * delta / chord, for the larger root
            const float root = std::abs(dot) * chordTimesMiss +
                               chordTimesD0Max * std::sqrt(dot * dot + delta);

            // [CUT] doublet curvature: the least any circle through the
            // doublet can have, delta / (chord * hypot(delta, root)), against
            // the loosest limit
            if ((doubletSlides && !curvatureSettled
                     ? curvatureCutSlidingSquare
                     : curvatureCutLoosestSquare) *
                    cChord * cChord * (delta * delta + root * root) <
                delta * delta) {
              continue;
            }
          }

          if (!tauSettled) {
            allTausSettled = false;
          }

          const float hypotTau = fastHypot(1, doubletTau);
          const float expEta = hypotTau - doubletTau;
          // 1 / expEta, since (hypotTau - tau) * (hypotTau + tau) == 1
          const float invExpEta = hypotTau + doubletTau;

          // [CUT] match before create: n2 must have an incoming edge that
          // agrees on tau, once it has enough edges to decide
          if (useMatchBeforeCreate) {
            bool isGood = n2NumEdges <= m_cfg.matchBeforeCreateMaxEdges;

            if (!isGood) {
              // n2's doublets are in exp(eta) order, so only a window of
              // them can match
              const auto [n2InFirst, expEtaHigh] =
                  tauRatioWindow(slw.doublets->expEta, n2FirstEdge, n2LastEdge,
                                 expEta, m_cfg.tauRatioPrecut);

              for (std::uint32_t n2InIdx = n2InFirst; n2InIdx < n2LastEdge;
                   ++n2InIdx) {
                // the doublet's own tau, as it has no fit of its own
                const float tau2 = slw.doublets->expEta[n2InIdx];

                if (tau2 > expEtaHigh) {
                  break;
                }

                const float tauRatio = tau2 * invExpEta - 1.0f;

                if (std::abs(tauRatio) > m_cfg.tauRatioPrecut) {  // bad match
                  continue;
                }
                isGood = true;  // good match found
                break;
              }
            }

            if (!isGood) {  // no match found, skip creating [n1 <- n2] edge
              continue;
            }
          }

          // [CUT] settled ends share a transverse position
          const float turnChord = fastHypot(turnDx, turnDy);
          if (turnChord == 0.0f) {
            continue;
          }

          // unit chord direction, for the turn cut and cached on the edge
          const float invChord = 1.0f / turnChord;
          const float cosAlpha12 = turnDx * invChord;
          const float sinAlpha12 = turnDy * invChord;

          // middle node of every triplet tried below
          const std::array<float, 3> point2{position2[0], position2[1], z2};

          // Kept whether or not a triplet passes through it now: it may yet
          // be the outer doublet of one made from a bin further in.
          inner.doublets.emplace_back(
              n2Idx, turnChord, cosAlpha12, sinAlpha12,
              static_cast<std::uint32_t>(inner.properties.size()),
              slw.connectionIdx, calibrate && slides2, turnSettled, tauSettled);
          inner.expEta.push_back(expEta);
          if (tauPull) {
            inner.etaGradient.push_back(
                settledCuts
                    ? chordEtaGradient(settled1, settled2, nodeStripAxis[n1Idx],
                                       nodeStripAxis[n2Idx])
                    : chordEtaGradient(point1, {position2[0], position2[1], z2},
                                       nodeStripAxis[n1Idx],
                                       nodeStripAxis[n2Idx]));
          }
          if (settledFit) {
            inner.settledEnds.push_back({settled1, settled2});
          }

          ++numCreatedEdges;
          ++nDoublets;

          // Both stay put for the loop below, which adds no doublet. Taken
          // after the one above went in, as a bin that links to itself has
          // `outer` and `inner` the same.
          detail::DisplacedGbtsDoublet& newDoublet = inner.doublets.back();
          const auto newDoubletIdx =
              static_cast<std::uint32_t>(inner.doublets.size() - 1);
          BinDoublets& outer = *slw.doublets;

          // n2's doublets are in exp(eta) order, so the loosest tau ratio
          // cut leaves a window of them, found by bisection
          // the tighter window holds only if no partner could be unsettled
          const bool windowSettled = settleTau && tauSettled && allTausSettled;
          const float windowTauRatioCut =
              slw.loosestTauRatioCut +
              (windowSettled ? slw.loosestTauRatioCorrStripSettled
                             : slw.loosestTauRatioCorrStrip);
          // With a pull, widen to the largest uncertainty any of n2's
          // doublets could give the pair; each candidate is held to its own
          // below.
          float windowRatio = windowTauRatioCut;
          std::array<float, 2> newGradient{};
          if (tauPull) {
            newGradient = inner.etaGradient.back();
            const std::array<float, 2>& most = maxEtaGradient[n2Idx];
            const float shared = std::abs(newGradient[1]) + most[0];
            const float bound2 = newGradient[0] * newGradient[0] +
                                 shared * shared + most[1] * most[1];
            const float pullRatio =
                std::expm1(m_cfg.tauRatioPullCut *
                           std::sqrt(stripRes2 * bound2 + pullFloor2));
            windowRatio =
                tauFixed ? std::max(windowTauRatioCut, pullRatio) : pullRatio;
          }
          const auto [inEdgeFirst, expEtaHigh] = tauRatioWindow(
              outer.expEta, n2FirstEdge, n2LastEdge, expEta, windowRatio);

          // The curvature limit the turn cut holds this doublet's triplets
          // to: the low-eta one where its tau keeps them all below the split.
          const bool lowEtaTurn =
              lowEtaTauBound > 0.0f && std::abs(doubletTau) < lowEtaTauBound;
          const float turnCurvatureCut =
              lowEtaTurn ? curvatureCutLowEta : curvatureCutLoosest;
          const float turnCurvatureCutSliding =
              lowEtaTurn ? curvatureCutLowEta * m_cfg.curvatureScaleEndcapStrip
                         : curvatureCutSliding;

          // looking for neighbours of the new doublet using outer node
          for (std::uint32_t inEdgeIdx = inEdgeFirst; inEdgeIdx < n2LastEdge;
               ++inEdgeIdx) {
            // [CUT] new doublet's triplet count full: an unrecorded triplet
            // could not be matched later anyway
            if (newDoublet.nProperties >=
                detail::kDisplacedGbtsMaxEdgeNeighbours) {
              break;
            }

            // [CUT] past the tau ratio window, and so is every later one
            if (outer.expEta[inEdgeIdx] > expEtaHigh) {
              break;
            }

            const float tauRatio = outer.expEta[inEdgeIdx] * invExpEta;
            const float absTauRatio = std::abs(tauRatio - 1.0f);

            // [CUT] loosest tau ratio, off the compact array. With a pull it
            // waits for the turn cut instead, which costs less than the pull's
            // square root and log and removes most of what reaches it.
            if (!tauPull && absTauRatio > windowTauRatioCut) {
              continue;
            }

            detail::DisplacedGbtsDoublet* pS = &outer.doublets[inEdgeIdx];

            // [CUT] outer edge's neighbour array full
            if (pS->nNei >= detail::kDisplacedGbtsMaxEdgeNeighbours) {
              continue;
            }

            // whether the third node slides radially, which loosens the cut
            const bool outerSlides = pS->n2Slides;

            const bool skipTurn =
                innerSkipTurn || !turnSettled || !pS->chordSettled;

            // Turn between the two chords, from the cached unit vectors.
            // For a circle sin(turn) = curvature * L13 <= curvature *
            // (L12 + L23). Per pair rather than hoisted: hoisting leaked 2%
            // more candidates to the fit for 4% more time.
            if (!skipTurn) {
              const float sinDelta =
                  cosAlpha12 * pS->sinAlpha - sinAlpha12 * pS->cosAlpha;
              const float cosDelta =
                  cosAlpha12 * pS->cosAlpha + sinAlpha12 * pS->sinAlpha;

              // [CUT] turn past a quarter circle, where the sine test below
              // stops meaning anything
              if (cosDelta <= 0.0f) {
                continue;
              }

              // [CUT] turn angle against the loosest curvature limit, looser
              // still where an end was settled on a radial strip
              if (std::abs(sinDelta) > ((doubletSlides || outerSlides)
                                            ? turnCurvatureCutSliding
                                            : turnCurvatureCut) *
                                           (turnChord + pS->chord)) {
                continue;
              }
            }

            // the pull's tolerance on the doublets' difference in eta
            float pullTolerance = 0.f;
            float logTauRatio = 0.f;
            if (tauPull) {
              pullTolerance =
                  m_cfg.tauRatioPullCut *
                  std::sqrt(stripRes2 *
                                tripletEtaGradient2(
                                    newGradient, outer.etaGradient[inEdgeIdx]) +
                            pullFloor2);
              logTauRatio = std::log(tauRatio);

              // [CUT] loosest tau ratio, with the pull
              if (!passesTau(absTauRatio, windowTauRatioCut, logTauRatio,
                             pullTolerance)) {
                continue;
              }
            }

            // the candidate triplet, inside out
            const std::array<SpacePointIndex, 3> tripletNodes{n1Idx, n2Idx,
                                                              pS->n2};

            // whether each node's tau carries an along-strip coordinate
            const std::array<bool, 3> isStrip{
                !isPixel1, !isPixel2,
                nodeView.strip(tripletNodes[2]) != nullptr};

            // the looser of the cuts the two layer pairs carry
            const GbtsLayerConnection& innerPair = *slw.connection;
            const GbtsLayerConnection& outerPair = connections[pS->connection];
            const float pairTauRatioCut = std::max(innerPair.tauRatioCutInner,
                                                   outerPair.tauRatioCutOuter);

            // Widen for strips, whose along-strip coordinate is unresolved
            // until the calibration below, which removes this again; less
            // where both doublets already settled theirs.
            const float stripTauRatioCorr =
                (isStrip[0] || isStrip[1] || isStrip[2])
                    ? ((settleTau && tauSettled && pS->tauSettled)
                           ? std::max(innerPair.tauRatioCorrStripSettled,
                                      outerPair.tauRatioCorrStripSettled)
                           : std::max(innerPair.tauRatioCorrStrip,
                                      outerPair.tauRatioCorrStrip))
                    : 0.0f;

            // [CUT] doublet tau ratio, loose; retaken tight after calibration
            if (!passesTau(absTauRatio, pairTauRatioCut + stripTauRatioCorr,
                           logTauRatio, pullTolerance)) {
              continue;
            }

            // The fit starts from the nominal nodes, or from where the two
            // doublets settled them: the inner end of the inner doublet, the
            // outer end of the outer one, and the middle node between the two
            // places its two chords cross its strip.
            std::array<std::array<float, 3>, 3> points{
                point1, point2, nodePoint(tripletNodes[2])};

            if (settledFit) {
              const auto& innerEnds = inner.settledEnds[newDoubletIdx];
              const auto& outerEnds = outer.settledEnds[inEdgeIdx];
              points[0] = innerEnds[0];
              for (std::uint32_t c = 0; c < 3; ++c) {
                points[1][c] = 0.5f * (innerEnds[1][c] + outerEnds[0][c]);
              }
              points[2] = outerEnds[1];
            }

            std::optional<detail::TripletCircle> circle =
                fitTripletCircle(points);

            // [CUT] degenerate circle fit
            if (!circle.has_value()) {
              continue;
            }

            // strip widening left over: kept unless every strip end resolved
            float unresolvedStripCorr = stripTauRatioCorr;

            // whether any end moved off its nominal position
            bool moved = false;

            // Slide each strip end to where the fitted track crosses it, then
            // refit once; the slide is small enough for one pass.
            if (calibrate) {
              bool resolved = true;
              bool allStripsResolved = true;

              for (std::uint32_t k = 0; k < 3; ++k) {
                const auto* strip = nodeView.strip(tripletNodes[k]);

                if (strip == nullptr) {
                  if (isStrip[k]) {
                    allStripsResolved = false;
                  }
                  continue;
                }

                if (!Acts::detail::calibrateOuterStripSpacePoint(
                        circle->direction(k), *strip, points[k],
                        m_cfg.maxStripLengthFraction)) {
                  resolved = false;
                  break;
                }

                moved = true;
              }

              // [CUT] fitted direction misses a strip
              if (!resolved) {
                continue;
              }

              if (moved) {
                circle = fitTripletCircle(points);

                // [CUT] degenerate refit
                if (!circle.has_value()) {
                  continue;
                }
              }

              if (allStripsResolved) {
                unresolvedStripCorr = 0.0f;
              }
            }

            // Tau ratio retaken on the calibrated ends at the tight
            // threshold. If nothing moved it equals the doublet one.
            float resolvedTauRatio = absTauRatio;
            float resolvedLogRatio = logTauRatio;
            float tightGradient2 = 0.f;

            if (tauPull) {
              tightGradient2 = tripletEtaGradient2(
                  newGradient, outer.etaGradient[inEdgeIdx]);
            }

            if (moved) {
              const std::optional<float> expEta12 =
                  chordExpEta(points[0], points[1]);
              const std::optional<float> expEta23 =
                  chordExpEta(points[1], points[2]);

              // [CUT] calibrated ends share a transverse position
              if (!expEta12.has_value() || !expEta23.has_value()) {
                continue;
              }

              resolvedTauRatio = std::abs(*expEta23 / *expEta12 - 1.0f);

              if (tauPull) {
                resolvedLogRatio = std::log(*expEta23 / *expEta12);
                const std::array<float, 3>& axis1 = nodeStripAxis[n1Idx];
                const std::array<float, 3>& axis2 = nodeStripAxis[n2Idx];
                const std::array<float, 3>& axis3 = nodeStripAxis[pS->n2];
                tightGradient2 = tripletEtaGradient2(
                    chordEtaGradient(points[0], points[1], axis1, axis2),
                    chordEtaGradient(points[1], points[2], axis2, axis3));
              }
            }

            const float tightPullTolerance =
                tauPull ? m_cfg.tauRatioPullCutTight *
                              std::sqrt(stripRes2 * tightGradient2 + pullFloor2)
                        : 0.f;

            // [CUT] tight tau ratio on the calibrated ends
            if (!passesTau(resolvedTauRatio,
                           pairTauRatioCut + unresolvedStripCorr,
                           resolvedLogRatio, tightPullTolerance)) {
              continue;
            }

            const float tripletCurv = circle->curvature;

            // [CUT] eta dependent curvature
            if (std::abs(tripletCurv) >
                (std::abs(circle->tau) < m_cfg.curvatureSplitAbsTau
                     ? curvatureCutLowEta
                     : curvatureCutHighEta)) {
              continue;
            }

            // With d0 fitted, z0 and z(rOut) are single values rather than
            // the doublet's bands, from the calibrated inner end.
            if (m_cfg.tripletFilterRZ) {
              const float r1c = fastHypot(points[0][0], points[0][1]);
              const std::optional<float> s1 =
                  arcFromPerigee(r1c, circle->d0, tripletCurv);

              // [CUT] triplet never reaches its own inner radius
              if (!s1.has_value()) {
                continue;
              }

              const float z0 = points[0][2] - *s1 * circle->tau;

              // [CUT] triplet z0
              if (z0 < m_cfg.minZ0 || z0 > m_cfg.maxZ0) {
                continue;
              }

              // a track whose |d0| exceeds rOut never gets there, and has
              // nothing to check
              if (const std::optional<float> sOut =
                      arcFromPerigee(rOut, circle->d0, tripletCurv);
                  sOut.has_value()) {
                const float zOuter = points[0][2] + (*sOut - *s1) * circle->tau;

                // [CUT] triplet z(rOut)
                if (zOuter < cutZMinU || zOuter > cutZMaxU) {
                  continue;
                }
              }
            }

            if (m_cfg.validateTriplets) {
              // [CUT] d0
              if (std::abs(circle->d0) > m_cfg.d0Max) {
                continue;
              }

              if (tripletCurv != 0.0f) {  // straight-line track is OK
                // curvature is 1 / 2R, so pT = bFieldInZ * R
                const float pT = 0.5f * std::abs(bFieldInZ / tripletCurv);

                // [CUT] min pT
                if (pT < tripletPtMin) {
                  continue;
                }

                // [CUT] tighter tau ratio for relatively high-pT tracks
                if (pT > 5 * tripletPtMin) {
                  if (!passesTau(resolvedTauRatio, 0.9f * m_cfg.tauRatioCut,
                                 resolvedLogRatio, 0.9f * tightPullTolerance)) {
                    continue;
                  }
                }
              }
            }

            // Match against the triplets the outer edge already belongs to,
            // which share the (n2, n3) doublet: the displaced form of the
            // prompt dphi/dcurv doublet match. An outer edge in no triplet
            // yet is a chain end and passes.
            if (pS->nProperties > 0) {
              // tangent at n2, where the recorded ones were taken
              const float phiShared = circle->tangentPhi(1);

              bool matched = false;

              for (std::uint8_t p = 0; p < pS->nProperties; ++p) {
                const detail::TripletProperties& prev =
                    outer.properties[pS->firstProperty + p];

                float dPhi = phiShared - prev.phi;

                if (dPhi < -std::numbers::pi_v<float>) {
                  dPhi += 2 * std::numbers::pi_v<float>;
                } else if (dPhi > std::numbers::pi_v<float>) {
                  dPhi -= 2 * std::numbers::pi_v<float>;
                }

                // [CUT] triplet match: tangent dphi
                if (std::abs(dPhi) > m_cfg.cutDPhiMax) {
                  continue;
                }

                const float dcurv = tripletCurv - prev.curvature;

                // [CUT] triplet match: dcurv
                if (dcurv < -m_cfg.cutDCurvMax || dcurv > m_cfg.cutDCurvMax) {
                  continue;
                }

                // [CUT] triplet match: tau ratio of the two fits
                if (std::abs(prev.expEta * circle->invExpEta - 1.0f) >
                    pairTauRatioCut + unresolvedStripCorr) {
                  continue;
                }

                matched = true;
                break;
              }

              // [CUT] no recorded triplet matched
              if (!matched) {
                continue;
              }
            }

            // The triplet passed, so both its doublets belong in the graph.
            const std::size_t nNewEdges =
                (pS->edge == detail::kNoGbtsEdge ? 1u : 0u) +
                (newDoublet.edge == detail::kNoGbtsEdge ? 1u : 0u);

            // [CUT] graph full
            if (edgeStorage.size() + nNewEdges > nMaxEdges) {
              graphFull = true;
              continue;
            }

            if (pS->edge == detail::kNoGbtsEdge) {
              pS->edge = static_cast<std::uint32_t>(edgeStorage.size());
              edgeStorage.emplace_back(n2Idx, pS->n2, outer.expEta[inEdgeIdx]);
            }

            if (newDoublet.edge == detail::kNoGbtsEdge) {
              newDoublet.edge = static_cast<std::uint32_t>(edgeStorage.size());
              edgeStorage.emplace_back(n1Idx, n2Idx, expEta);
            }

            // record the tangent at n1 for later triplets through (n1, n2)
            inner.properties.emplace_back(circle->expEta, tripletCurv,
                                          circle->tangentPhi(0));
            ++newDoublet.nProperties;

            if (newDoublet.nProperties ==
                detail::kDisplacedGbtsMaxEdgeNeighbours) {
              ++nSaturatedProperties;
            }

            // the outer edge takes the inner one as a neighbour
            detail::DisplacedGbtsEdge& outerEdge = edgeStorage[pS->edge];

            outerEdge.vNei[outerEdge.nNei] = newDoublet.edge;
            ++outerEdge.nNei;
            ++pS->nNei;

            if (pS->nNei == detail::kDisplacedGbtsMaxEdgeNeighbours) {
              ++nSaturatedNeighbours;
            }

            nConnections++;
          }  // inEdgeIdx
        }  // n2PhiIdx
      }  // slw

      // Put n1's doublets in exp(eta) order, so that whatever matches
      // against them reads only the window its tau ratio cut leaves. Each
      // doublet carries where its triplet properties and its edge are, and
      // nothing else holds a doublet by its place, so moving them is safe.
      if (numCreatedEdges > 1) {
        const std::uint32_t first = edgeInfo[n1Idx].firstEdge;

        sortOrder.resize(numCreatedEdges);
        for (std::uint32_t k = 0; k < numCreatedEdges; ++k) {
          sortOrder[k] = first + k;
        }
        std::ranges::sort(sortOrder, {}, [&inner](const std::uint32_t idx) {
          return inner.expEta[idx];
        });

        sortedDoublets.clear();
        sortedExpEta.clear();
        sortedEtaGradient.clear();
        sortedSettledEnds.clear();
        for (const std::uint32_t idx : sortOrder) {
          sortedDoublets.push_back(inner.doublets[idx]);
          sortedExpEta.push_back(inner.expEta[idx]);
          if (tauPull) {
            sortedEtaGradient.push_back(inner.etaGradient[idx]);
          }
          if (settledFit) {
            sortedSettledEnds.push_back(inner.settledEnds[idx]);
          }
        }
        std::ranges::copy(sortedDoublets, inner.doublets.begin() + first);
        std::ranges::copy(sortedExpEta, inner.expEta.begin() + first);
        if (tauPull) {
          std::ranges::copy(sortedEtaGradient,
                            inner.etaGradient.begin() + first);
        }
        if (settledFit) {
          std::ranges::copy(sortedSettledEnds,
                            inner.settledEnds.begin() + first);
        }
      }

      // the largest rates of n1's doublets, for the windows of the doublets
      // that will end on n1
      if (tauPull && numCreatedEdges > 0) {
        const std::uint32_t first = edgeInfo[n1Idx].firstEdge;
        std::array<float, 2>& most = maxEtaGradient[n1Idx];
        for (std::uint32_t k = first; k < first + numCreatedEdges; ++k) {
          most[0] = std::max(most[0], std::abs(inner.etaGradient[k][0]));
          most[1] = std::max(most[1], std::abs(inner.etaGradient[k][1]));
        }
      }

      // makes n1's doublets visible to the next bin group in
      edgeInfo[n1Idx].numEdges = numCreatedEdges;
    }  // n1Idx

    finishGroup(bg);
  }  // bin groups

  const auto nEdges = static_cast<std::uint32_t>(edgeStorage.size());

  ACTS_DEBUG("Made " << nDoublets << " doublets, holding at most "
                     << maxHeldDoublets << " at once, and kept " << nEdges
                     << " as graph edges");

  if (graphFull) {
    ACTS_WARNING(
        "Maximum number of graph edges exceeded - possible efficiency loss "
        << nEdges);
  }

  if (nSaturatedNeighbours > 0) {
    ACTS_WARNING("Maximum number of edge connections ("
                 << detail::kDisplacedGbtsMaxEdgeNeighbours << ") reached on "
                 << nSaturatedNeighbours << " of " << nEdges
                 << " edges - possible efficiency loss");
  }

  if (nSaturatedProperties > 0) {
    ACTS_WARNING("Maximum number of recorded triplets ("
                 << detail::kDisplacedGbtsMaxEdgeNeighbours << ") reached on "
                 << nSaturatedProperties << " of " << nEdges
                 << " edges - possible efficiency loss");
  }

  graph.nDoublets = nDoublets;
  graph.nConnections = nConnections;
  return graph;
}

}  // namespace Acts::Experimental
