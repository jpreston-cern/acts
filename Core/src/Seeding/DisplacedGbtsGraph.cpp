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

namespace Acts::Experimental{

namespace {

/// Sliding window in phi used to define range used for edge creation.
///
/// Covers one non-empty source eta bin, whose phi-ordered node list it holds
/// directly so that the innermost loop does not reach through the bin.
struct SlidingWindow {
  /// phi-ordered nodes of the bin, including the wrap-around duplicates
  const std::pair<float, SpacePointIndex>* phiNodes{};
  /// number of entries in @c phiNodes
  std::uint32_t numPhiNodes{};
  /// sliding window position
  std::uint32_t firstIt{};
  /// window half-width;
  float deltaPhi{};
  /// Inside-out pixel barrel ordinal of the bin's layer, -1 for the rest.
  std::int32_t barrelOrder{-1};
  /// Type of the bin's layer.
  GbtsLayerType type{};
  /// Technology of the bin's layer.
  GbtsLayerTechnology technology{};
};

/// exp(eta) of the chord between two points, the form the tau ratio of a
/// triplet compares its two doublets in.
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

}  // namespace
    
  DisplacedGbtsGraph::DisplacedGbtsGraph(const Config& config,
                     std::shared_ptr<const GbtsGeometry> geometry,
                     std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(config),
      m_geometry(std::move(geometry)),
      m_logger(std::move(logger)) {}

  std::pair<std::uint32_t, std::uint32_t> DisplacedGbtsGraph::buildTheGraph(
    const GbtsRoiDescriptor& roi, GbtsNodeStorage& nodeStorage,
    std::vector<detail::DisplacedGbtsEdge>& edgeStorage, const float bFieldInZ) const {
  // used to calculate the outer Z cut on doublets
  const float cutZMinU =
      m_cfg.minZ0 + m_cfg.maxOuterRadius * static_cast<float>(roi.dzdrMin());
  const float cutZMaxU =
      m_cfg.maxZ0 + m_cfg.maxOuterRadius * static_cast<float>(roi.dzdrMax());

  // these define the actual pT used 
  const float tripletPtMin = m_cfg.tripletPtFraction * m_cfg.minPt;

  const float ptScale = m_cfg.tuningPt / m_cfg.minPt;
  
  // eta dependent curvature cuts, on triplets now that a doublet no longer
  // knows its own curvature
  // different values due to different affects of MS at differnet pTs
  const float curvatureCutHighEta = m_cfg.maxCurvatureHighEta * ptScale;
  const float curvatureCutLowEta = m_cfg.maxCurvatureLowEta * ptScale;

  // the looser of the two, for the precuts that run before the fitted tau
  // exists to say which of them applies
  const float curvatureCutLoosest =
      std::max(curvatureCutLowEta, curvatureCutHighEta);


  // the loosest tau ratio threshold the triplet matching can apply
  const float maxTauRatioCut =
      m_cfg.tauRatioCut + (m_cfg.useAdaptiveCuts ? m_cfg.tauRatioCorr : 0.0f) +
      (nodeStorage.hasStrips() ? m_cfg.tauRatioCorrStrip : 0.0f);

  // the default sliding window along phi. Taken from the node storage so that
  // the windows and the phi indexing they slide over cannot disagree.
  const float deltaPhi0 = 0.5f * nodeStorage.m_cfg.phiSliceWidth;

  std::uint32_t nConnections = 0;

  // edges that filled their neighbour array, each one a connection the graph
  // could have gone on to make and now cannot
  std::uint32_t nSaturatedNeighbours = 0;

  // edges that filled their triplet property store, which leaves them unable
  // to record what a later edge would have been matched against
  std::uint32_t nSaturatedProperties = 0;

  edgeStorage.reserve(m_cfg.nMaxEdges);
  
  // number of edges acepted into the storage, 
  // will need a seperate one for the proper storage, 
  // global and local layer storages should be seperate due to some properties not being needed 
  std::uint32_t nEdges = 0;

  // views of the nodes and edges (need to doublec check the update on the edge view, this should be keyed to global, not local layer pair store)
  const detail::GbtsNodeView nodeView = nodeStorage.nodeView();
  const std::span<const detail::GbtsNodeParams> params =
      nodeStorage.nodeParams();
  const std::span<detail::GbtsNodeEdgeInfo> edgeInfo =
      nodeStorage.nodeEdgeInfo();

  // reused across bin groups so that the windows are allocated once
  std::vector<SlidingWindow> phiSlidingWindow;
  
  // if we apply calibration to the strips
  const bool calibrate = m_cfg.calibrateStrips && nodeStorage.hasStrips();

  // the nominal position of a node, which is where a strip node sits before
  // its triplet resolves it along the strip
  const auto nodePoint = [&nodeView](const SpacePointIndex node) {
    const std::array<float, 4>& position = nodeView.positions[node];
    return std::array<float, 3>{position[0], position[1], position[2]};
  };

  // loop over bin groups
  for (const auto& bg : m_geometry->binGroups()) {
    // B1 is the innermost bin in the pair, 
    // but we still look out to in on the bin group level,
    // has to be this way to set the width of the phi window (deifned from innermost bin)
    const detail::GbtsEtaBinInfo& B1 = nodeStorage.etaBin(bg.bin);

    if (B1.empty()) {
      continue;
    }

    // used for the phi window width creation
    const float rb1 = B1.minRadius;
    // still aplicacable but will be for strips this time
    const std::int32_t barrelOrder1 = B1.barrelOrder;

    // used to define whether a node in this bin needs calibrating due to low stip resoution in the r/ plane
    const bool isPixel1 = B1.technology == GbtsLayerTechnology::Pixel;
    // defines whether we are a pixel barrel layer, used for triplet validation (not needed),
    // and the tau correction (which should be used)
    const bool isPixelBarrel1 = barrelOrder1 >= 0;
    
    // this should be kept, but we shouldnt use "isConnected" and "zbitmask"
    // due to them not making much sense on a min level = 2 graph
    const bool useMatchBeforeCreate =
        m_cfg.matchBeforeCreate && barrelOrder1 >= 0 &&
        barrelOrder1 <= m_cfg.matchBeforeCreateMaxBarrelOrder;

    // prepare a sliding window for each non-empty bin2 in the group

    phiSlidingWindow.clear();

    // ready for loop over outer bins 
    // loop over bin groups
  // loop over n2 eta-bins in L2 layers
    for (const std::uint32_t b2Idx : bg.links) {
      // outer most bin 
      const detail::GbtsEtaBinInfo& B2 = nodeStorage.etaBin(b2Idx);

      if (B2.empty()) {
        continue;
      }

      // used for working out phi window width
      const float rb2 = B2.maxRadius;

      float deltaPhi = deltaPhi0;  // the default

      // override the default window width
      if (m_cfg.useEtaBinning) {
        
        const float absDr = std::fabs(rb2 - rb1);
        const float maxD0 = m_cfg.d0Max;
        auto phiWindow = [&rb1, &rb2, &maxD0, ptScale](const float& phiWindowOffset, const float& phiWindowSlope){
                                   const float maxD0Square = maxD0*maxD0;

                                   // How far round the beamline a displaced
                                   // track walks between the two radii. A
                                   // straight track of impact parameter d0
                                   // sits at azimuth acos(d0 / r), which is
                                   // exact and not a small angle expansion.
                                   //
                                   // A track only reaches radius r if its |d0|
                                   // is below r, so a bin inside the d0 limit
                                   // takes the limit down to its own radius.
                                   // Without that the ratio passes one and the
                                   // arc cosine returns a NaN, which would
                                   // leave every phi comparison false and so
                                   // the sliding window neither closing nor
                                   // advancing.
                                   const float frac1 = std::min(1.0f, maxD0/rb1);
                                   const float frac2 = std::min(1.0f, maxD0/rb2);
                                   const float displacmentTerm = std::acos(frac2) - std::acos(frac1);

                                   const float corr1 = std::sqrt(std::max(0.0f, rb1*rb1 - maxD0Square));
                                   const float corr2 = std::sqrt(std::max(0.0f, rb2*rb2 - maxD0Square));

                                   const float curvatureTerm = (corr2 - corr1)*phiWindowSlope*ptScale;

                                   // The two add. Where the track started and
                                   // which way it bends are independent, so a
                                   // track displaced to one side and bending
                                   // the same way walks the sum of them, and a
                                   // window that took the difference would be
                                   // covering only the case where they cancel.
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
      // again this is problamatic ( think i need to implement a layer first approach here)
      // maybe keep for now and let claude do it once you have 
      window.barrelOrder = B2.barrelOrder;
      window.type = B2.type;
      window.technology = B2.technology;
    }

    // in GBTSv3 the outer loop goes over n1 nodes in the Layer 1 bin (inner most node first)
    for (SpacePointIndex n1Idx = B1.nodes.first; n1Idx < B1.nodes.second;
         ++n1Idx) {
      
      // put inside loop here where things are going to get more complicated 
      // initialization using the top watermark of the edge storage
      edgeInfo[n1Idx].firstEdge = nEdges;

      // the counter for the incoming graph edges created for n1
      // it has to be here as we iterrate over edges on a shared node in the inner loop
      std::uint16_t numCreatedEdges = 0;
      
      const detail::GbtsNodeParams& n1pars = params[n1Idx];

      const float phi1 = n1pars.phi;
      const float r1 = n1pars.r;
      const float z1 = n1pars.z;

      // the chord of a pair, both paths need it now
      const std::array<float, 4>& position = nodeView.positions[n1Idx];
      float x1 = position[0];
      float y1 = position[1];
      
      // the intermediate loop over sliding windows. these are associated with the nodes in the outer bin
      for (auto& slw : phiSlidingWindow) {
        const std::int32_t barrelOrder2 = slw.barrelOrder;

        const bool isPixel2 = slw.technology == GbtsLayerTechnology::Pixel;
        const bool isPixelBarrel2 = barrelOrder2 >= 0;

        const float deltaPhi = slw.deltaPhi;

        // sliding window phi1 +/- deltaPhi
        // need to check if the equation i derived is the half width or full phi window
        const float minPhi = phi1 - deltaPhi;
        const float maxPhi = phi1 + deltaPhi;

        for (std::uint32_t n2PhiIdx = slw.firstIt; n2PhiIdx < slw.numPhiNodes;
             ++n2PhiIdx) {

          const float phi2 = slw.phiNodes[n2PhiIdx].first;

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
          

          // On the nominal radii, so an endcap pair the slide would have
          // opened up is lost here. The cheap reject is worth more.

          // we should work out the chord here i think for each pair and then cut on a minimum chord length 
          // set to same value i think for now as they serve similar purposes (can tune later)
          if (chord < m_cfg.minDeltaRadius) {
            continue;
          }

          const float z2 = n2pars.z;

          // The ends stay nominal here. Resolving a strip node needs the
          // track direction, and the chord is only that for a track from the
          // beamline, so the slide along the strip waits for the triplet fit
          // below. What that leaves behind is a doublet tau taken on
          // unresolved ends, which is what `tauRatioCorrStrip` is for.
          const float dz = z2 - z1;
          const float tau = dz / chord;
          const float ftau = std::fabs(tau);
          if (ftau > m_cfg.maxAbsTau) {
            continue;
          }

          // z is linear in transverse arc length, so z0 = z1 - S1 * tau with S1
          // the transverse path from the perigee to hit 1. The prompt form took
          // S1 = r1, which only holds for a track leaving the beamline
          // radially. Displaced, S1 = sqrt(r1^2 - d0^2), and d0 is unknown at
          // doublet stage -- two transverse points do not fix a circle -- so z0
          // is a band spanned by d0 in [0, d0Max]. S1 is monotonic in d0, so
          // the two endpoints bracket the band.
          const float d0MaxSquare = m_cfg.d0Max * m_cfg.d0Max;

          const float s1Max = r1;  // d0 = 0
          const float s1Min = std::sqrt(std::max(0.f, r1 * r1 - d0MaxSquare));

          const float z0A = z1 - s1Max * tau;
          const float z0B = z1 - s1Min * tau;
          const float z0Lo = std::min(z0A, z0B);
          const float z0Hi = std::max(z0A, z0B);

          if (m_cfg.doubletFilterRZ) {
            // the band has to overlap the allowed range, not sit inside it
            if (z0Hi < m_cfg.minZ0 || z0Lo > m_cfg.maxZ0) {
              continue;
            }

            // z(rOut) = z1 + (sOut - S1) * tau, and (sOut - S1) also grows
            // monotonically with d0, so the endpoints bracket this band too.
            const float rOut = m_cfg.maxOuterRadius;
            const float sOutMin =
                std::sqrt(std::max(0.f, rOut * rOut - d0MaxSquare));

            const float zOutA = z1 + (rOut - s1Max) * tau;
            const float zOutB = z1 + (sOutMin - s1Min) * tau;

            // NOTE: the two bands are tested separately although they share the
            // same d0, so a pair can pass the two on different d0 values. That
            // is a relaxation, never a rejection, which is the safe direction
            // for a pre-filter; the triplet fit resolves d0 properly.
            if (std::max(zOutA, zOutB) < cutZMinU ||
                std::min(zOutA, zOutB) > cutZMaxU) {
              continue;
            }
          }

        
          const float hypotTau = fastHypot(1, tau);
          const float expEta = hypotTau - tau;
          // 1 / expEta, since (hypotTau - tau) * (hypotTau + tau) == 1. The
          // sum is also the better conditioned form for large tau.
          const float invExpEta = hypotTau + tau;

          // match edge candidate against edges incoming to n2
          // can keep this i think 
          if (useMatchBeforeCreate) {
            // we must have enough incoming edges to decide
            bool isGood = n2NumEdges <= m_cfg.matchBeforeCreateMaxEdges;

            if (!isGood) {
              const float uat1 = invExpEta;

              for (std::uint32_t n2InIdx = n2FirstEdge; n2InIdx < n2LastEdge;
                   ++n2InIdx) {
                // the doublet's own tau: a displaced edge has no fit
                // parameters of its own until a triplet gives it some
                const float tau2 = edgeStorage[n2InIdx].expEta;
                const float tauRatio = tau2 * uat1 - 1.0f;

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

          // The inner edge's own chord direction, which every candidate below
          // is measured against and which its own edge then carries for the
          // triplets it will be the outer edge of. The chord is already here,
          // so the unit vector costs a reciprocal.
          const float invChord = 1.0f / chord;
          const float cosAlpha12 = dx * invChord;
          const float sinAlpha12 = dy * invChord;

          if (nEdges < m_cfg.nMaxEdges) {
            edgeStorage.emplace_back(n1Idx, n2Idx, expEta, chord, cosAlpha12,
                                     sinAlpha12, barrelOrder2);

            ++numCreatedEdges;

            const std::uint32_t outEdgeIdx = nEdges;
            
            const float uat2 = invExpEta;
          
            // looking for neighbours of the new edge using outer node 
            for (std::uint32_t inEdgeIdx = n2FirstEdge; inEdgeIdx < n2LastEdge;
                 ++inEdgeIdx) {

              // The new edge can only record so many triplets, and a triplet
              // it did not record is one a later match cannot see, so there is
              // nothing to gain from looking further.
              if (edgeStorage[outEdgeIdx].properties.size() >=
                  detail::kGbtsMaxEdgeNeighbours) {
                break;
              }

              // i think i need to add the triplet stuff here already such that everything can be calibrated properly before (maybe)
              detail::DisplacedGbtsEdge* pS = &edgeStorage[inEdgeIdx];

              const float absTauRatio = std::abs(pS->expEta * uat2 - 1.0f);

              // rejects most candidates before the layer bookkeeping below
              // this can stay as we want the tau ratio still 
              if (absTauRatio > maxTauRatioCut) {
                continue;
              }

              if (pS->nNei >= detail::kGbtsMaxEdgeNeighbours) {
                continue;
              }

              // The three nodes of the candidate triplet, inside out. The
              // doublet stage could not fix a circle, so everything the prompt
              // graph read off a single doublet -- curvature, the tangent
              // azimuth, d0, pT -- is found here instead.
              const std::array<SpacePointIndex, 3> tripletNodes{n1Idx, n2Idx,
                                                                pS->n2};

              // Whether each node has an along-strip coordinate sitting in the
              // tau it contributes. The inner two go by their bin, the third
              // by whether it carries a stereo pair.
              const std::array<bool, 3> isStrip{
                  !isPixel1, !isPixel2,
                  nodeView.strip(tripletNodes[2]) != nullptr};

              // A strip end has still to slide along its strip, and in the
              // endcap that slide is largely radial, so it moves the chords
              // the two geometric cuts below are reading. Where they cannot be
              // sure they do not cut: rejecting a candidate the fit would have
              // kept is the one thing they must not do.
              const bool endsWillMove =
                  calibrate && (isStrip[0] || isStrip[1] || isStrip[2]);

              // The turn from the inner edge to this one, which is
              // asin(curvature * L13) and so falls with pT. Both directions
              // were taken once, each when its own edge was made, so the turn
              // comes out of a dot and a cross of two cached unit vectors and
              // the cut never forms the angle at all.
              //
              // L13 is not measured here; the two chords bound it, since
              // L13 <= L12 + L23. That bound is this pair's own rather than
              // the detector's, which is the difference between a cut that
              // fires and one that does not: the widest turn the curvature
              // limit allows runs from about three degrees between adjacent
              // pixel layers to ten across the strips.
              if (!endsWillMove) {
                const float sinDelta =
                    cosAlpha12 * pS->sinAlpha - sinAlpha12 * pS->cosAlpha;
                const float cosDelta =
                    cosAlpha12 * pS->cosAlpha + sinAlpha12 * pS->sinAlpha;

                // Turning by a quarter circle between two layers puts the
                // triplet's own ends half a circle apart, which is orders
                // below any pT worth seeding. It also puts the turn past where
                // its sine still grows with curvature, so the test below would
                // stop meaning anything.
                if (cosDelta <= 0.0f) {
                  continue;
                }

                if (std::abs(sinDelta) >
                    curvatureCutLoosest * (chord + pS->chord)) {
                  continue;
                }
              }

              const std::int32_t barrelOrder3 = pS->n2BarrelOrder;

              const bool isPixelBarrel3 = barrelOrder3 >= 0;

              float addTauRatioCorr = 0;

              if (m_cfg.useAdaptiveCuts) {
                if (isPixelBarrel1 && isPixelBarrel2 && isPixelBarrel3) {
                  // three radially consecutive layers, none skipped
                  const bool noGap = (barrelOrder2 - barrelOrder1) == 1 &&
                                     (barrelOrder3 - barrelOrder2) == 1;

                  // assume more scattering due to the layer in between
                  if (!noGap) {
                    addTauRatioCorr = m_cfg.tauRatioCorr;
                  }
                } else {
                  bool mixedTriplet =
                      isPixelBarrel1 && isPixelBarrel2 && !isPixelBarrel3;
                  if (mixedTriplet) {
                    addTauRatioCorr = m_cfg.tauRatioCorr;
                  }
                }
              }
              // The two doublets sharing a strip node resolved it separately,
              // so a triplet through a strip may disagree on tau by more. Any
              // of the three: the outer two carry their end's error into tau.
              // Held apart from the correction above because the calibration
              // below is what takes this one away and not that one.
              const float stripTauRatioCorr =
                  (m_cfg.tauRatioCorrStrip > 0.f &&
                   (isStrip[0] || isStrip[1] || isStrip[2]))
                      ? m_cfg.tauRatioCorrStrip
                      : 0.0f;

              // bad match. Loose on purpose: the along-strip coordinate the
              // doublet tau leans on is the very thing the calibration is
              // about to move, so this is only a precut and the tight version
              // of it runs on the resolved points further down.
              if (absTauRatio >
                  m_cfg.tauRatioCut + addTauRatioCorr + stripTauRatioCorr) {
                continue;
              }

              std::array<std::array<float, 3>, 3> points{
                  nodePoint(tripletNodes[0]), nodePoint(tripletNodes[1]),
                  nodePoint(tripletNodes[2])};

              std::optional<detail::TripletCircle> circle = fitTripletCircle(points);

              if (!circle.has_value()) {
                continue;
              }

              // How much of the strip looseness survives. It comes off only
              // the nodes the calibration actually put back where the track
              // crossed, so a node its bin calls a strip but that carries no
              // stereo pair keeps it, as does every node when the calibration
              // is off.
              float unresolvedStripCorr = stripTauRatioCorr;

              // whether any end ended up somewhere other than where the node
              // storage put it, which is the only case the cuts below have
              // anything new to work with
              bool moved = false;

              // Slide each strip node along its strip to where the fitted
              // track crossed it, then fit again on the moved ends. Once is
              // enough: the slide is a small fraction of a strip and the fit
              // is linear in the node positions to that order.
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

                  // a direction that misses the strip is no crossing at all
                  if (!Acts::detail::calibrateOuterStripSpacePoint(
                          circle->direction(k), *strip, points[k],
                          m_cfg.maxStripLengthFraction)) {
                    resolved = false;
                    break;
                  }

                  moved = true;
                }

                if (!resolved) {
                  continue;
                }

                if (moved) {
                  circle = fitTripletCircle(points);

                  if (!circle.has_value()) {
                    continue;
                  }
                }

                if (allStripsResolved) {
                  unresolvedStripCorr = 0.0f;
                }
              }

              // The along-strip coordinate is where a strip node's z lives and
              // z is all tau is made of, which is why the doublet tau above
              // needed the looseness at all. Every strip end now sits where
              // the fitted track crossed it, so the same ratio is worth
              // retaking on the resolved points, at the tight threshold. It
              // still says something despite the fit having supplied the
              // directions: an end moves by about the module separation per
              // unit of tau error, a millimetre against a lever arm of a
              // hundred, so the calibration cannot pull three ends into an
              // agreement a real track would not have had.
              //
              // Nothing moved means this is the doublet tau over again, and
              // the threshold is the one it already passed, so it is only
              // retaken when it can differ.
              float resolvedTauRatio = absTauRatio;

              if (moved) {
                const std::optional<float> expEta12 =
                    chordExpEta(points[0], points[1]);
                const std::optional<float> expEta23 =
                    chordExpEta(points[1], points[2]);

                if (!expEta12.has_value() || !expEta23.has_value()) {
                  continue;
                }

                resolvedTauRatio = std::abs(*expEta23 / *expEta12 - 1.0f);
              }

              if (resolvedTauRatio >
                  m_cfg.tauRatioCut + addTauRatioCorr + unresolvedStripCorr) {
                continue;
              }

              const float tripletCurv = circle->curvature;

              // the eta dependent curvature cut the doublets could not carry
              if (std::abs(tripletCurv) >
                  (std::abs(circle->tau) < m_cfg.curvatureSplitAbsTau
                       ? curvatureCutLowEta
                       : curvatureCutHighEta)) {
                continue;
              }

              // final check: cuts on pT and d0
              if (m_cfg.validateTriplets) {
                if (std::abs(circle->d0) > m_cfg.d0Max) {
                  continue;
                }

                if (tripletCurv != 0.0f) {  // straight-line track is OK
                  // curvature is 1 / R in the prompt convention, where R is
                  // twice the radius, so pT = bFieldInZ * R / 2 is this
                  const float pT = 0.5f * std::abs(bFieldInZ / tripletCurv);

                  if (pT < tripletPtMin) {
                    continue;
                  }

                  if (pT > 5 * tripletPtMin) {  // relatively high-pT track

                    if (resolvedTauRatio > 0.9f * m_cfg.tauRatioCut) {
                      continue;
                    }
                  }
                }
              }

              // Match against the triplets the outer edge is already part of.
              // Each of those shares the (n2, n3) doublet with this one, so
              // the two are the same track candidate seen one node apart and
              // both hold a tangent azimuth at n2: this is the displaced
              // stand-in for the dphi/dcurv match the prompt graph runs
              // between two doublets, which needed a curvature a displaced
              // doublet does not have. An outer edge that is not yet in any
              // triplet is the outermost pair of a chain and has nothing to
              // disagree with, so it passes.
              if (!pS->properties.empty()) {
                // node 1 of this triplet is the inner node of the shared
                // doublet, which is where the recorded tangents were taken
                const float phiShared = circle->tangentPhi(1);

                bool matched = false;

                for (const detail::TripletProperties& prev : pS->properties) {
                  float dPhi = phiShared - prev.phi;

                  if (dPhi < -std::numbers::pi_v<float>) {
                    dPhi += 2 * std::numbers::pi_v<float>;
                  } else if (dPhi > std::numbers::pi_v<float>) {
                    dPhi -= 2 * std::numbers::pi_v<float>;
                  }

                  if (std::abs(dPhi) > m_cfg.cutDPhiMax) {
                    continue;
                  }

                  const float dcurv = tripletCurv - prev.curvature;

                  if (dcurv < -m_cfg.cutDCurvMax || dcurv > m_cfg.cutDCurvMax) {
                    continue;
                  }

                  // the arc corrected taus of the two fits, which is a tighter
                  // statement than the doublet tau ratio already tested
                  if (std::abs(prev.expEta * circle->invExpEta - 1.0f) >
                      m_cfg.tauRatioCut + addTauRatioCorr +
                          unresolvedStripCorr) {
                    continue;
                  }

                  matched = true;
                  break;
                }

                if (!matched) {
                  continue;
                }
              }

              // The new edge is the inner edge of this triplet, so what a
              // later triplet through (n1, n2) will be matched against is the
              // tangent at n1, node 0.
              edgeStorage[outEdgeIdx].properties.emplace_back(
                  circle->expEta, tripletCurv, circle->tangentPhi(0));

              if (edgeStorage[outEdgeIdx].properties.size() ==
                  detail::kGbtsMaxEdgeNeighbours) {
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

      // updating the n1 node attributes. Without this the edges just created
      // are invisible to the nodes of the next bin group in and no triplet
      // could ever be formed through n1.
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

  return std::make_pair(nEdges, nConnections);
}

std::uint32_t DisplacedGbtsGraph::runCCA(
    const std::uint32_t nEdges,
    std::vector<detail::DisplacedGbtsEdge>& edgeStorage) const {
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
    std::vector<detail::DisplacedGbtsEdge>& edgeStorage,
    std::uint32_t nEdges) const {
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
