#include "Acts/Seeding/DisplacedGbtsGraph.hpp"
#include "Acts/SpacePointFormation/detail/StripSpacePointCalibrationImpl.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
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
  
  // eta dependent curvature cuts on doublets
  // different values due to different affects of MS at differnet pTs
  const float curvatureCutHighEta = m_cfg.maxCurvatureHighEta * ptScale;
  const float curvatureCutLowEta = m_cfg.maxCurvatureLowEta * ptScale;

  // the loosest tau ratio threshold the triplet matching can apply
  const float maxTauRatioCut =
      m_cfg.tauRatioCut + (m_cfg.useAdaptiveCuts ? m_cfg.tauRatioCorr : 0.0f) +
      (nodeStorage.hasStrips() ? m_cfg.tauRatioCorrStrip : 0.0f);

  // the default sliding window along phi. Taken from the node storage so that
  // the windows and the phi indexing they slide over cannot disagree.
  const float deltaPhi0 = 0.5f * nodeStorage.m_cfg.phiSliceWidth;

  std::uint32_t nConnections = 0;

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

  // Put a strip node back where a direction says it crossed. Both ends, since
  // the nominal position sits on the inner strip and the calibrated one on the
  // outer. A direction that misses either strip is no crossing at all.
  const auto calibrateNode = [&](const SpacePointIndex node,
                                 const std::array<float, 3>& direction,
                                 float& r, float& z) {
    std::array<float, 3> point{};
    if (!Acts::detail::calibrateOuterStripSpacePoint(
            direction, nodeStorage.strip(node), point,
            m_cfg.maxStripLengthFraction)) {
      return false;
    }
    r = fastHypot(point[0], point[1]);
    z = point[2];
    return true;
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
        // the general equations for this need to change
        auto phiWindow = [&rb1, &rb2, &maxD0, ptScale](const float& phiWindowOffset, const float& phiWindowSlope){
                                   const float maxD0Square = maxD0*maxD0;
                                   const float rb1Square = rb1*rb1;
                                   const float rb2Square = rb2*rb2;

                                   const float frac1 = maxD0/rb1;
                                   const float frac2 = maxD0/rb2;
                                   const float displacmentTerm = std::acos(frac2) - std::acos(frac1);

                                   const float corr1 = std::sqrt(rb1Square - maxD0Square);
                                   const float corr2 = std::sqrt(rb2Square - maxD0Square); 
                                   
                                   const float curvatureTerm = (corr2 - corr1)*phiWindowSlope*ptScale;

                                   const float absPhiWindow = phiWindowOffset + std::abs(displacmentTerm - curvatureTerm);

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

        const bool stripPair = calibrate && (!isPixel1 || !isPixel2);

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

          // get rid of for now, 
          // even though we are only using one large storage atm, we should move to two in the future
          // which wont require this
          const std::uint16_t nodeInfo = n2Info.isConnected;

          const std::uint32_t n2FirstEdge = n2Info.firstEdge;
          const std::uint16_t n2NumEdges = n2Info.numEdges;
          const std::uint32_t n2LastEdge = n2FirstEdge + n2NumEdges;

          const detail::GbtsNodeParams& n2pars = params[n2Idx];

          const std::array<float, 4>& position = nodeView.positions[n2Idx];

          const float dx = position[0] - x1;
          const float dy = position[1] - y1;

          const float r2 = n2pars.r;
          const float chord = fastHypot(dx, dy);
          

          // On the nominal radii, so an endcap pair the slide would have
          // opened up is lost here. The cheap reject is worth more.

          // we should work out the chord here i think for each pair and then cut on a minimum chord length 
          // set to same value i think for now as they serve similar purposes (can tune later)
          if (chord < m_cfg.minDeltaRadius) {
            continue;
          }

          const float z2 = n2pars.z;

          // the ends as the pair puts them: nominal, or slid along the strip
          // when resolved. Azimuth is kept as it was -- exactly right in the
          // barrel, and to the stereo angle in the endcap, where the strip
          // slid along is a hair off radial.
          float r1c = r1;
          float z1c = z1;
          float r2c = r2;
          float z2c = z2;

          const float dz = z2c - z1c;
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

          const float s1Max = r1c;  // d0 = 0
          const float s1Min = std::sqrt(std::max(0.f, r1c * r1c - d0MaxSquare));

          const float z0A = z1c - s1Max * tau;
          const float z0B = z1c - s1Min * tau;
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

            const float zOutA = z1c + (rOut - s1Max) * tau;
            const float zOutB = z1c + (sOutMin - s1Min) * tau;

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
                const float tau2 = edgeStorage[n2InIdx].properties[0].expEta;
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

          if (nEdges < m_cfg.nMaxEdges) {
            edgeStorage.emplace_back(n1Idx, n2Idx, expEta, barrelOrder2);

            ++numCreatedEdges;

            const std::uint32_t outEdgeIdx = nEdges;
            
            const float uat2 = invExpEta;
          
            // FIXME: incomplete. In GbtsGraph this block also runs the triplet
            // neighbour loop over n2's incoming edges (the tau/dPhi/dcurv
            // matching that fills vNei and sets isConnected) and ends with
            // ++nEdges. Until that is written, uat2/phi2u/curv2 are unused and
            // nEdges never advances, so every edge overwrites slot 0.
            // looking for neighbours of the new edge using outer node 
            for (std::uint32_t inEdgeIdx = n2FirstEdge; inEdgeIdx < n2LastEdge;
                 ++inEdgeIdx) {
              
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
              if (m_cfg.tauRatioCorrStrip > 0.f &&
                  (!isPixel1 || !isPixel2 ||
                   nodeView.strip(pS->n2) != nullptr)) {
                addTauRatioCorr += m_cfg.tauRatioCorrStrip;
              }

              // bad match
              if (absTauRatio > m_cfg.tauRatioCut + addTauRatioCorr) {
                continue;
              }
            } // inEdgeIdx
          } // if (nEdges < nMaxEdges)
        } // n2PhiIdx
      } // slw
    } // n1Idx
  } // bin groups

  return std::make_pair(nEdges, nConnections);
}

} // Acts::Experimental namespace