#include "Acts/Seeding/DisplacedGbtsGraph.hpp"

#include "Acts/SpacePointFormation/detail/StripSpacePointCalibrationImpl.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <array>
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
    std::vector<detail::GbtsEdge>& edgeStorage, const float bFieldInZ) const {
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
    
    // this should be kept, but we shouldnt "isConnected" and "zbitmask"
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
        // the general equations for this need to change
        if (absDr < m_cfg.phiWindowSplitDeltaRadius) {
          deltaPhi = m_cfg.phiWindowNearOffset +
                     m_cfg.phiWindowNearSlope * ptScale * absDr;
        } else {
          deltaPhi = m_cfg.phiWindowFarOffset +
                     m_cfg.phiWindowFarSlope * ptScale * absDr;
        }
      }

      SlidingWindow& window = phiSlidingWindow.emplace_back();
      window.phiNodes = B2.phiNodes.data();
      window.numPhiNodes = static_cast<std::uint32_t>(B2.phiNodes.size());
      window.deltaPhi = deltaPhi;
      // again this is problamatic ( think i need to implement a layer first approach here)
      window.barrelOrder = B2.barrelOrder;
      window.type = B2.type;
      window.technology = B2.technology;
    }

    // in GBTSv3 the outer loop goes over n1 nodes in the Layer 1 bin
    for (SpacePointIndex n1Idx = B1.nodes.first; n1Idx < B1.nodes.second;
         ++n1Idx) {

      // put inside loop here where things are going to get more complicated 
    }
  }

  return std::make_pair(nEdges, nConnections);
}

} // Acts::Experimental namespace