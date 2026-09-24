#include "Acts/Seeding/DisplacedGbtsGraph.hpp"
#include "Acts/SpacePointFormation/detail/StripSpacePointCalibrationImpl.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <numbers>
#include <optional>
#include <span>
#include <utility>
#include <vector>

// Every rejection in buildTheGraph is marked `[CUT]`, so `grep "\[CUT\]"`
// lists them in the order a candidate meets them.

namespace Acts::Experimental{

namespace {

/// Sliding window in phi over one non-empty eta bin, holding its phi-ordered
/// nodes directly so the innermost loop does not reach through the bin.
struct SlidingWindow {
  /// phi-ordered nodes of the bin, including the wrap-around duplicates
  const std::pair<float, SpacePointIndex>* phiNodes{};
  /// number of entries in @c phiNodes
  std::uint32_t numPhiNodes{};
  /// sliding window position
  std::uint32_t firstIt{};
  /// window half-width
  float deltaPhi{};
  /// How deep the bin's layer sits in the barrel, -1 for an endcap.
  std::int32_t depth{-1};
  /// Technology of the bin's layer.
  GbtsLayerTechnology technology{};
  /// Type of the bin's layer.
  GbtsLayerType type{};
};

/// Whether a strip node on this layer slides radially when calibrated: true
/// for endcap strips. Barrel strips slide along z and keep their transverse
/// position.
constexpr bool slidesRadially(GbtsLayerType type,
                              GbtsLayerTechnology technology) {
  return type == GbtsLayerType::Endcap &&
         technology == GbtsLayerTechnology::Strip;
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

}  // namespace

  DisplacedGbtsGraph::DisplacedGbtsGraph(const Config& config,
                     std::shared_ptr<const GbtsGeometry> geometry,
                     std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(config),
      m_geometry(std::move(geometry)),
      m_logger(std::move(logger)) {}

  detail::GbtsGraph<detail::DisplacedGbtsEdge> DisplacedGbtsGraph::buildTheGraph(
    const GbtsRoiDescriptor& roi, GbtsNodeStorage& nodeStorage,
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

  // the loosest tau ratio threshold the triplet matching can apply
  const float maxTauRatioCut =
      m_cfg.tauRatioCut + (m_cfg.useAdaptiveCuts ? m_cfg.tauRatioCorr : 0.0f) +
      (nodeStorage.hasStrips() ? m_cfg.tauRatioCorrStrip : 0.0f);

  // per-event parts of the doublet z0 band over d0 in [0, d0Max]
  const float d0MaxSquare = m_cfg.d0Max * m_cfg.d0Max;
  const float rOut = m_cfg.maxOuterRadius;
  const float sOutMin = std::sqrt(std::max(0.f, rOut * rOut - d0MaxSquare));

  // default phi half-window, from the node storage so the two cannot disagree
  const float deltaPhi0 = 0.5f * nodeStorage.m_cfg.phiSliceWidth;

  std::uint32_t nConnections = 0;

  // edges that filled their neighbour array or their triplet property store,
  // reported as possible efficiency loss
  std::uint32_t nSaturatedNeighbours = 0;
  std::uint32_t nSaturatedProperties = 0;

  edgeStorage.reserve(m_cfg.nMaxEdges);

  // Every edge's exp(eta), kept beside the storage so the first tau ratio cut
  // scans 4 bytes per edge instead of loading whole edges. Same indices.
  std::vector<float> edgeExpEta;
  edgeExpEta.reserve(m_cfg.nMaxEdges);

  std::uint32_t nEdges = 0;

  const detail::GbtsNodeView nodeView = nodeStorage.nodeView();
  const std::span<const detail::GbtsNodeParams> params =
      nodeStorage.nodeParams();
  const std::span<detail::GbtsNodeEdgeInfo> edgeInfo =
      nodeStorage.nodeEdgeInfo();

  // reused across bin groups so that the windows are allocated once
  std::vector<SlidingWindow> phiSlidingWindow;

  // if we apply calibration to the strips
  const bool calibrate = m_cfg.calibrateStrips && nodeStorage.hasStrips();

  // nominal position of a node, before any strip calibration
  const auto nodePoint = [&nodeView](const SpacePointIndex node) {
    const std::array<float, 4>& position = nodeView.positions[node];
    return std::array<float, 3>{position[0], position[1], position[2]};
  };

  // loop over bin groups
  for (const auto& bg : m_geometry->binGroups()) {
    // B1 is the inner bin of every pair in the group
    const detail::GbtsEtaBinInfo& B1 = nodeStorage.etaBin(bg.bin);

    if (B1.empty()) {
      continue;
    }

    // used for the phi window width
    const float rb1 = B1.minRadius;

    const bool isPixel1 = B1.technology == GbtsLayerTechnology::Pixel;
    // barrel depth over all barrel layers, -1 for an endcap
    const std::int32_t depth1 = B1.depth;

    const bool useMatchBeforeCreate = m_cfg.matchBeforeCreate && depth1 >= 0 &&
                                      depth1 <= m_cfg.matchBeforeCreateMaxDepth;

    // prepare a sliding window for each non-empty bin2 in the group
    phiSlidingWindow.clear();

    for (const std::uint32_t b2Idx : bg.links) {
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
        auto phiWindow = [&rb1, &rb2, &maxD0, ptScale](const float& phiWindowOffset, const float& phiWindowSlope){
                                   const float maxD0Square = maxD0*maxD0;

                                   // Azimuth a straight track of impact
                                   // parameter d0 walks between the radii,
                                   // exactly acos(d0 / r). The ratio is capped
                                   // at 1 for bins inside d0Max, otherwise
                                   // acos gives NaN and the window breaks.
                                   const float frac1 = std::min(1.0f, maxD0/rb1);
                                   const float frac2 = std::min(1.0f, maxD0/rb2);
                                   const float displacmentTerm = std::acos(frac2) - std::acos(frac1);

                                   const float corr1 = std::sqrt(std::max(0.0f, rb1*rb1 - maxD0Square));
                                   const float corr2 = std::sqrt(std::max(0.0f, rb2*rb2 - maxD0Square));

                                   const float curvatureTerm = (corr2 - corr1)*phiWindowSlope*ptScale;

                                   // displacement and bending are independent,
                                   // so the worst case is their sum
                                   const float absPhiWindow = phiWindowOffset + std::abs(displacmentTerm) + std::abs(curvatureTerm);

                                   return absPhiWindow;
                                   };
        if (absDr < m_cfg.phiWindowSplitDeltaRadius) {
          deltaPhi = phiWindow(m_cfg.phiWindowNearOffset, m_cfg.phiWindowNearSlope);
        } else {
          deltaPhi = phiWindow(m_cfg.phiWindowFarOffset, m_cfg.phiWindowFarSlope);
        }
      }

      SlidingWindow& window = phiSlidingWindow.emplace_back();
      window.phiNodes = B2.phiNodes.data();
      window.numPhiNodes = static_cast<std::uint32_t>(B2.phiNodes.size());
      window.deltaPhi = deltaPhi;
      window.depth = B2.depth;
      window.technology = B2.technology;
      window.type = B2.type;
    }

    // outer loop over n1 nodes in the inner bin
    for (SpacePointIndex n1Idx = B1.nodes.first; n1Idx < B1.nodes.second;
         ++n1Idx) {
      // initialization using the top watermark of the edge storage
      edgeInfo[n1Idx].firstEdge = nEdges;

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

      // Transverse path from the perigee to n1 at d0 = d0Max, and the lever
      // arms out to rOut at both ends of the d0 range, for the z0 band.
      const float s1Min = std::sqrt(std::max(0.f, r1 * r1 - d0MaxSquare));
      const float outerLeverZeroD0 = rOut - r1;
      const float outerLeverMaxD0 = sOutMin - s1Min;

      // loop over the sliding windows of the outer bins
      for (auto& slw : phiSlidingWindow) {
        const std::int32_t depth2 = slw.depth;

        const bool isPixel2 = slw.technology == GbtsLayerTechnology::Pixel;

        // The turn cut is skipped when switched off, or when an end will
        // slide radially under calibration (an endcap strip). Settled here for
        // the inner two nodes; the third is checked per candidate.
        const bool innerSkipTurn =
            !m_cfg.useTurnAngleCut ||
            (calibrate && (slidesRadially(B1.type, B1.technology) ||
                           slidesRadially(slw.type, slw.technology)));

        // deltaPhi is the half width: displacement and charge are both
        // symmetric about phi1
        const float deltaPhi = slw.deltaPhi;
        const float minPhi = phi1 - deltaPhi;
        const float maxPhi = phi1 + deltaPhi;

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

          const SpacePointIndex n2Idx = slw.phiNodes[n2PhiIdx].second;

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
          // graph's radial step cut; nominal positions.
          if (chord < m_cfg.minDeltaRadius) {
            continue;
          }

          const float z2 = n2pars.z;

          // Doublet tau on nominal ends: strip calibration needs the track
          // direction, which only the triplet fit gives. `tauRatioCorrStrip`
          // covers the difference.
          const float dz = z2 - z1;
          const float tau = dz / chord;
          const float ftau = std::fabs(tau);

          // [CUT] max |tau|
          if (ftau > m_cfg.maxAbsTau) {
            continue;
          }

          // [CUT] cluster width tau window, on both nodes. Only pixel barrel
          // nodes have a finite window; the rest pass.
          if (ftau < n1pars.minTau || ftau > n1pars.maxTau) {
            continue;
          }

          if (ftau < n2pars.minTau || ftau > n2pars.maxTau) {
            continue;
          }

          // z0 = z1 - S1 * tau, with S1 = sqrt(r1^2 - d0^2) the transverse path
          // from the perigee. d0 is unknown for a doublet, so z0 and z(rOut)
          // are bands over d0 in [0, d0Max], bracketed by the two ends.
          if (m_cfg.doubletFilterRZ) {
            const float z0A = z1 - r1 * tau;     // d0 = 0
            const float z0B = z1 - s1Min * tau;  // d0 = d0Max

            // [CUT] doublet z0 band must overlap [minZ0, maxZ0]
            if (std::max(z0A, z0B) < m_cfg.minZ0 ||
                std::min(z0A, z0B) > m_cfg.maxZ0) {
              continue;
            }

            const float zOutA = z1 + outerLeverZeroD0 * tau;
            const float zOutB = z1 + outerLeverMaxD0 * tau;

            // [CUT] doublet z(rOut) band must overlap the RoI. The two bands
            // may pass on different d0; `tripletFilterRZ` tightens this.
            if (std::max(zOutA, zOutB) < cutZMinU ||
                std::min(zOutA, zOutB) > cutZMaxU) {
              continue;
            }
          }

          const float hypotTau = fastHypot(1, tau);
          const float expEta = hypotTau - tau;
          // 1 / expEta, since (hypotTau - tau) * (hypotTau + tau) == 1
          const float invExpEta = hypotTau + tau;

          // [CUT] match before create: n2 must have an incoming edge that
          // agrees on tau, once it has enough edges to decide
          if (useMatchBeforeCreate) {
            bool isGood = n2NumEdges <= m_cfg.matchBeforeCreateMaxEdges;

            if (!isGood) {
              for (std::uint32_t n2InIdx = n2FirstEdge; n2InIdx < n2LastEdge;
                   ++n2InIdx) {
                // the doublet's own tau, as the edge has no fit of its own
                const float tau2 = edgeStorage[n2InIdx].expEta;
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

          // unit chord direction, for the turn cut and cached on the edge
          const float invChord = 1.0f / chord;
          const float cosAlpha12 = dx * invChord;
          const float sinAlpha12 = dy * invChord;

          // middle node of every triplet tried below
          const std::array<float, 3> point2{position2[0], position2[1], z2};

          // [CUT] edge storage full
          if (nEdges < m_cfg.nMaxEdges) {
            edgeStorage.emplace_back(n1Idx, n2Idx, expEta, chord, cosAlpha12,
                                     sinAlpha12, depth2, slw.type,
                                     slw.technology);
            edgeExpEta.push_back(expEta);

            ++numCreatedEdges;

            const std::uint32_t outEdgeIdx = nEdges;

            // safe to hold: the storage was reserved to nMaxEdges
            detail::DisplacedGbtsEdge& newEdge = edgeStorage[outEdgeIdx];

            // looking for neighbours of the new edge using outer node
            for (std::uint32_t inEdgeIdx = n2FirstEdge; inEdgeIdx < n2LastEdge;
                 ++inEdgeIdx) {

              // [CUT] new edge's triplet store full: an unrecorded triplet
              // could not be matched later anyway
              if (newEdge.nProperties >= detail::kGbtsMaxEdgeNeighbours) {
                break;
              }

              const float absTauRatio =
                  std::abs(edgeExpEta[inEdgeIdx] * invExpEta - 1.0f);

              // [CUT] loosest tau ratio, off the compact array
              if (absTauRatio > maxTauRatioCut) {
                continue;
              }

              detail::DisplacedGbtsEdge* pS = &edgeStorage[inEdgeIdx];

              // [CUT] outer edge's neighbour array full
              if (pS->nNei >= detail::kGbtsMaxEdgeNeighbours) {
                continue;
              }

              const bool skipTurn =
                  innerSkipTurn ||
                  (calibrate && slidesRadially(pS->n2Type, pS->n2Technology));

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

                // [CUT] turn angle against the loosest curvature limit
                if (std::abs(sinDelta) >
                    curvatureCutLoosest * (chord + pS->chord)) {
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

              const std::int32_t depth3 = pS->n2Depth;

              float addTauRatioCorr = 0;

              // Widen for material: a gap in barrel depth (over all barrel
              // layers, strips included), or a third node that left the barrel.
              if (m_cfg.useAdaptiveCuts && depth1 >= 0 && depth2 >= 0) {
                if (depth3 >= 0) {
                  const bool noGap =
                      (depth2 - depth1) == 1 && (depth3 - depth2) == 1;

                  if (!noGap) {
                    addTauRatioCorr = m_cfg.tauRatioCorr;
                  }
                } else {
                  addTauRatioCorr = m_cfg.tauRatioCorr;
                }
              }

              // Widen for strips, whose along-strip coordinate is unresolved
              // until the calibration below, which removes this again.
              const float stripTauRatioCorr =
                  (m_cfg.tauRatioCorrStrip > 0.f &&
                   (isStrip[0] || isStrip[1] || isStrip[2]))
                      ? m_cfg.tauRatioCorrStrip
                      : 0.0f;

              // [CUT] doublet tau ratio, loose; retaken tight after calibration
              if (absTauRatio >
                  m_cfg.tauRatioCut + addTauRatioCorr + stripTauRatioCorr) {
                continue;
              }

              std::array<std::array<float, 3>, 3> points{
                  point1, point2, nodePoint(tripletNodes[2])};

              std::optional<detail::TripletCircle> circle = fitTripletCircle(points);

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
              }

              // [CUT] tight tau ratio on the calibrated ends
              if (resolvedTauRatio >
                  m_cfg.tauRatioCut + addTauRatioCorr + unresolvedStripCorr) {
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
                    if (resolvedTauRatio > 0.9f * m_cfg.tauRatioCut) {
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
                  const detail::TripletProperties& prev = pS->properties[p];

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
                      m_cfg.tauRatioCut + addTauRatioCorr +
                          unresolvedStripCorr) {
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

              // record the tangent at n1 for later triplets through (n1, n2)
              newEdge.properties[newEdge.nProperties] =
                  detail::TripletProperties{circle->expEta, tripletCurv,
                                            circle->tangentPhi(0)};
              ++newEdge.nProperties;

              if (newEdge.nProperties == detail::kGbtsMaxEdgeNeighbours) {
                ++nSaturatedProperties;
              }

              pS->vNei[pS->nNei] = outEdgeIdx;
              ++pS->nNei;

              if (pS->nNei == detail::kGbtsMaxEdgeNeighbours) {
                ++nSaturatedNeighbours;
              }

              nConnections++;
            } // inEdgeIdx
            nEdges++;
          } // if (nEdges < nMaxEdges)
        } // n2PhiIdx
      } // slw

      // makes n1's edges visible to the next bin group in
      edgeInfo[n1Idx].numEdges = numCreatedEdges;
    } // n1Idx
  } // bin groups

  if (nEdges >= m_cfg.nMaxEdges) {
    ACTS_WARNING(
        "Maximum number of graph edges exceeded - possible efficiency loss "
        << nEdges);
  }

  if (nSaturatedNeighbours > 0) {
    ACTS_WARNING("Maximum number of edge connections ("
                 << detail::kGbtsMaxEdgeNeighbours << ") reached on "
                 << nSaturatedNeighbours << " of " << nEdges
                 << " edges - possible efficiency loss");
  }

  if (nSaturatedProperties > 0) {
    ACTS_WARNING("Maximum number of recorded triplets ("
                 << detail::kGbtsMaxEdgeNeighbours << ") reached on "
                 << nSaturatedProperties << " of " << nEdges
                 << " edges - possible efficiency loss");
  }

  graph.nEdges = nEdges;
  graph.nConnections = nConnections;
  return graph;
}

std::uint32_t DisplacedGbtsGraph::runCCA(
    detail::GbtsGraph<detail::DisplacedGbtsEdge>& graph) const {
  const std::uint32_t nEdges = graph.nEdges;
  std::vector<detail::DisplacedGbtsEdge>& edgeStorage = graph.edgeStorage;

  std::uint32_t maxLevel = 0;

  std::uint32_t iter = 0;

  std::vector<detail::DisplacedGbtsEdge*> vOld;

  for (std::uint32_t edgeIndex = 0; edgeIndex < nEdges; ++edgeIndex) {
    detail::DisplacedGbtsEdge* pS = &(edgeStorage[edgeIndex]);
    if (pS->nNei == 0) {
      continue;
    }

    vOld.push_back(pS);
  }

  std::vector<detail::DisplacedGbtsEdge*> vNew;
  vNew.reserve(vOld.size());

  // generate proposals
  for (; iter < m_cfg.ccaMaxIterations; iter++) {
    vNew.clear();

    for (detail::DisplacedGbtsEdge* pS : vOld) {
      std::int32_t nextLevel = pS->level;

      for (std::uint32_t nIdx = 0; nIdx < pS->nNei; ++nIdx) {
        const std::uint32_t nextEdgeIdx = pS->vNei[nIdx];

        const detail::DisplacedGbtsEdge* pN = &(edgeStorage[nextEdgeIdx]);

        if (pS->level == pN->level) {
          nextLevel = pS->level + 1;
          vNew.push_back(pS);
          break;
        }
      }

      // proposal
      pS->next = static_cast<std::int8_t>(nextLevel);
    }

    // update

    std::uint32_t nChanges = 0;

    for (auto pS : vNew) {
      if (pS->next != pS->level) {
        nChanges++;
        pS->level = pS->next;
        // levels only grow from zero here, so the cast is safe
        maxLevel = std::max(maxLevel, static_cast<std::uint32_t>(pS->level));
      }
    }

    if (nChanges == 0) {
      break;
    }

    vOld.swap(vNew);
    vNew.clear();
  }

  return maxLevel;
}

std::vector<detail::DisplacedGbtsEdge*> DisplacedGbtsGraph::extractChainHeads(
    detail::GbtsGraph<detail::DisplacedGbtsEdge>& graph) const {
  const std::uint32_t nEdges = graph.nEdges;
  std::vector<detail::DisplacedGbtsEdge>& edgeStorage = graph.edgeStorage;

  const auto minLevel = static_cast<std::uint8_t>(m_cfg.minSeedLevel);
  // `addTriplets` accepts a chain one level short. Signed: an uncollected
  // edge sits at level -1 and `minSeedLevel` may be configured to 0.
  const int minLevelAddTriplets = int{minLevel} - 1;
  std::vector<detail::DisplacedGbtsEdge*> vChainHeads;

  vChainHeads.reserve(nEdges / 2);

  for (std::uint32_t edgeIndex = 0; edgeIndex < nEdges; ++edgeIndex) {
    detail::DisplacedGbtsEdge* pS = &edgeStorage[edgeIndex];

    if (!m_cfg.addTriplets) {
      if (pS->level < minLevel) {
        continue;
      }
    } else {  // eta-dependent cut
      // the edge's own eta, from the doublet: eta is -log(exp(eta)) and the
      // displaced edge keeps exp(eta) where the prompt one kept it in p[0]
      const float edgeAbsEta = std::abs(-std::log(pS->expEta));

      if (edgeAbsEta > m_cfg.maxAbsEtaAddTriplets) {
        if (pS->level < minLevel) {
          continue;
        }
      } else {
        if (pS->level < minLevelAddTriplets) {
          continue;
        }
      }
    }

    vChainHeads.push_back(pS);
  }

  if (vChainHeads.empty()) {
    return vChainHeads;
  }

  std::ranges::sort(vChainHeads, std::ranges::greater{},
                    [](const detail::DisplacedGbtsEdge* e) { return e->level; });

  return vChainHeads;
}

std::optional<detail::TripletCircle> DisplacedGbtsGraph::fitTripletCircle(
    const std::array<std::array<float, 3>, 3>& points) const {
  detail::TripletCircle circle{};

  const float x0 = points[1][0];
  const float y0 = points[1][1];
  const float r0 = fastHypot(x0, y0);

  if (r0 == 0.0f) {
    return std::nullopt;
  }

  const float cosA = x0 / r0;
  const float sinA = y0 / r0;

  circle.phiMid = std::atan2(y0, x0);

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

} // Acts::Experimental namespace
