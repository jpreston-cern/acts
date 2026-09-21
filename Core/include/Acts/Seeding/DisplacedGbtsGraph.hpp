
#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Seeding/GbtsNodeStorage.hpp"
#include "Acts/Seeding/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding/detail/DisplacedGraphTypes.hpp"
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

      /// Cut a triplet on the pT and d0 of its own fit. The fit itself runs
      /// either way, since the curvature properties an edge carries only
      /// exist once three nodes have fixed the circle.
      bool validateTriplets = true;

      /// Widens allowed variation in tau ratio if a layer is missed in edge
      /// connecting.
      bool useAdaptiveCuts = true;

      /// Tau ratio cut threshold.
      float tauRatioCut = 0.007f;

      /// Tau ratio precut threshold.
      float tauRatioPrecut = 0.009f;

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

      /// Apply RZ cuts on doublets.
      bool doubletFilterRZ = true;

      /// Maximum number of GBTS edges/doublets.
      std::uint32_t nMaxEdges = 2000000;

      /// Minimum delta radius between layers.
      float minDeltaRadius = 2.0f * Acts::UnitConstants::mm;

      /// Largest |cot(theta)| accepted for a doublet. The default corresponds to
      /// |eta| of about 4.3, beyond the acceptance of any current tracker.
      float maxAbsTau = 36.0f;

      /// Minimum z0 value. In pixel mode the value is picked from the RoI.
      float minZ0 = -600.0f;

      /// Maximum z0 value. In pixel mode the value is picked from the RoI.
      float maxZ0 = 600.0f;

      /// Maximum d0 impact parameter when validating an edge-connection triplet.
      float d0Max = 3.0f * Acts::UnitConstants::mm;

      /// Maximum difference in allowed tangent between candidate edge
      /// connections. Displaced, the two tangents being compared are the ones
      /// the two triplet fits put at the shared doublet's inner node.
      float cutDPhiMax = 0.012f;

      /// Maximum allowed curvature tolerance for candidate edge connections.
      float cutDCurvMax = 0.001f;

      /// pT at which the default cut coefficients were tuned; they scale by
      /// `tuningPt / minPt`.
      float tuningPt = 0.9f * Acts::UnitConstants::GeV;

      /// Maximum |curvature| above `curvatureSplitAbsTau`, before that scaling.
      float maxCurvatureHighEta = 4.75e-4f / Acts::UnitConstants::mm;

      /// Maximum |curvature| below `curvatureSplitAbsTau`, before that scaling.
      float maxCurvatureLowEta = 3.75e-4f / Acts::UnitConstants::mm;

      /// |cot(theta)| separating the two curvature cuts, corresponding to |eta|
      /// of about 2.1.
      float curvatureSplitAbsTau = 4.0f;

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

      /// Incoming edge count below which a node is accepted without a tau match.
      std::uint32_t matchBeforeCreateMaxEdges = 2;

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

      /// Maximum number of connected-component iterations.
      std::uint32_t ccaMaxIterations = 15;

      // Chain selection options, shared with the seed extraction that reads the
      // chains back out of the graph.

      /// Chain length a seed candidate must reach. Two, where the prompt graph
      /// asks for three: a link here already cost a triplet fit, so two edges
      /// are three nodes that were fitted and cut on rather than three that
      /// merely agreed on tau.
      std::uint32_t minSeedLevel = 2;

      /// optionally add 3 sp seeds within a certain eta range
      ///
      /// @note This accepts a chain one level short, so at the `minSeedLevel`
      ///       of two above it reaches down to a single edge: two nodes, no
      ///       fit, nothing this graph stands for. Raise `minSeedLevel` to
      ///       three before turning it on.
      bool addTriplets = false;

      /// the maximum allowed eta value in which
      /// three spacepoint seeds are passed through
      float maxAbsEtaAddTriplets = 1.5;
    };

    // we should give the same kind of objects in, just different configurations
    DisplacedGbtsGraph(const Config& config, std::shared_ptr<const GbtsGeometry> geometry,
                       std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
                           "DisplacedGbtsGraph", Acts::Logging::Level::INFO));

    const Config& config() const { return m_cfg; };

    std::pair<std::uint32_t, std::uint32_t> buildTheGraph(const GbtsRoiDescriptor& roi, 
                                                          GbtsNodeStorage& nodeStorage,
                                                          std::vector<detail::DisplacedGbtsEdge>& edgeStorage, 
                                                          float bFieldInZ) const;
            
    /// Run connected component analysis on the graph.
    /// @param nEdges Number of edges in the graph
    /// @param edgeStorage Storage containing graph edges
    /// @return The highest chain level any edge reached
    std::uint32_t runCCA(
        std::uint32_t nEdges,
        std::vector<detail::DisplacedGbtsEdge>& edgeStorage) const;

    /// extract edges that start a chain
    /// @param edgeStorage Storage containing graph edges
    /// @param nEdges Number of edges in the graph
    /// @return The edges that start chains, ordered by length of chain, empty
    ///         if no chain reached the required level
    std::vector<detail::DisplacedGbtsEdge*> extractChainHeads(
        std::vector<detail::DisplacedGbtsEdge>& edgeStorage,
        std::uint32_t nEdges) const;

    private:

    Config m_cfg;

    std::shared_ptr<const GbtsGeometry> m_geometry;

    std::unique_ptr<const Acts::Logger> m_logger;

    const Acts::Logger& logger() const { return *m_logger; }
  };

  
}