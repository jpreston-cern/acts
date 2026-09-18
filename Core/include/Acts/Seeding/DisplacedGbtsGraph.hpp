
#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Seeding/GbtsNodeStorage.hpp"
#include "Acts/Seeding/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace Acts::Experimental{

  class DisplacedGbtsGraph{

    public:

    struct Config{
      // going to fill as i go with this

      /// Match seeds before creating them.
      bool matchBeforeCreate = false;

      /// Widens allowed variation in tau ratio if a layer is missed in edge
      /// connecting.
      bool useAdaptiveCuts = true;

      /// Tau ratio cut threshold.
      float tauRatioCut = 0.007f;

      /// Correction applied to tau acceptance if a layer is missed during edge
      /// connecting.
      float tauRatioCorr = 0.006f;

      /// The same for a triplet any of whose three nodes a strip module made,
      /// whose two doublets resolved the shared node's along-strip coordinate
      /// separately. Reaches nothing without a strip in the triplet.
      float tauRatioCorrStrip = 0.03f;

      /// Minimum transverse momentum.
      float minPt = 1.0f * Acts::UnitConstants::GeV;

      /// Fraction of `minPt` a triplet may fall to, allowing for three-point
      /// pT resolution.
      float tripletPtFraction = 0.8f;

      // Graph-building options

      /// Use eta binning from geometry structure.
      bool useEtaBinning = true;

      /// Maximum number of GBTS edges/doublets.
      std::uint32_t nMaxEdges = 2000000;

      /// Minimum z0 value. In pixel mode the value is picked from the RoI.
      float minZ0 = -600.0f;

      /// Maximum z0 value. In pixel mode the value is picked from the RoI.
      float maxZ0 = 600.0f;

      /// Maximum d0 impact parameter when validating an edge-connection triplet.
      float d0Max = 3.0f * Acts::UnitConstants::mm;

      /// pT at which the default cut coefficients were tuned; they scale by
      /// `tuningPt / minPt`.
      float tuningPt = 0.9f * Acts::UnitConstants::GeV;

      /// Maximum |curvature| above `curvatureSplitAbsTau`, before that scaling.
      float maxCurvatureHighEta = 4.75e-4f / Acts::UnitConstants::mm;

      /// Maximum |curvature| below `curvatureSplitAbsTau`, before that scaling.
      float maxCurvatureLowEta = 3.75e-4f / Acts::UnitConstants::mm;

      /// Radial separation splitting the two phi-window slopes below.
      float phiWindowSplitDeltaRadius = 60.0f * Acts::UnitConstants::mm;

      /// Phi window below `phiWindowSplitDeltaRadius`, as an offset plus a slope
      /// times the radial separation.
      float phiWindowNearOffset = 0.002f;

      /// Slope of the near phi window per unit radial separation. Scaled by
      /// `tuningPt / minPt`.
      float phiWindowNearSlope = 4.33e-4f / Acts::UnitConstants::mm;

      /// Phi window above `phiWindowSplitDeltaRadius`, in the same form.
      float phiWindowFarOffset = 0.015f;

      /// Slope of the far phi window per unit radial separation. Scaled by
      /// `tuningPt / minPt`.
      float phiWindowFarSlope = 2.2e-4f / Acts::UnitConstants::mm;

      /// Highest pixel barrel layer, counted inside out, whose nodes are cut
      /// against the z0 histogram of their outer neighbourhood and whose isolated
      /// nodes are skipped. A negative value disables the cut.
      std::int32_t z0HistogramMaxBarrelOrder = 0;

      /// Highest pixel barrel layer, counted inside out, to which
      /// `matchBeforeCreate` applies when it is enabled. A negative value
      /// disables it.
      std::int32_t matchBeforeCreateMaxBarrelOrder = 1;

      /// Maximum radius of the pixel detector.
      float maxOuterRadius = 550.0f;

      /// Resolve a doublet's strip endpoints against its own direction before
      /// cutting on them. Nothing is written back; the correction belongs to the
      /// pair.
      bool calibrateStrips = true;

      /// How far along a strip a crossing may land and still be recovered, as a
      /// multiple of the strip half-length, so 1 is the strip itself. This is the
      /// same quantity as `TripletSeedFinder::Config::toleranceParam`.
      float maxStripLengthFraction = 1.1f;
    };

    // we should give the same kind of objects in, just different configurations
    DisplacedGbtsGraph(const Config& config, std::shared_ptr<const GbtsGeometry> geometry,
                       std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
                           "DisplacedGbtsGraph", Acts::Logging::Level::INFO));

    const Config& config() const { return m_cfg; };

    std::pair<std::uint32_t, std::uint32_t> buildTheGraph(const GbtsRoiDescriptor& roi, 
                                                          GbtsNodeStorage& nodeStorage,
                                                          std::vector<detail::GbtsEdge>& edgeStorage, 
                                                          float bFieldInZ) const;
            
    std::uint32_t runCCA(std::uint32_t nEdges, std::vector<detail::GbtsEdge>& edgeStorage) const;
    
    std::vector<detail::GbtsEdge*> extractChainHeads(std::vector<detail::GbtsEdge>& edgeStorage, std::uint32_t nEdges) const;

    private:

    Config m_cfg;

    std::shared_ptr<const GbtsGeometry> m_geometry;

    std::unique_ptr<const Acts::Logger> m_logger;

    const Acts::Logger& logger() const { return *m_logger; }
  };

  
}