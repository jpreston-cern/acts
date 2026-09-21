
#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsGraphConfig.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Seeding/GbtsNodeStorage.hpp"
#include "Acts/Seeding/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding/detail/DisplacedGraphTypes.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace Acts::Experimental{

  /// The doublet graph of the GBTS workflow, built without assuming the track
  /// came from the beamline.
  ///
  /// Does the same job as @c GbtsGraph and takes the same configuration, but
  /// a doublet here cannot be read for a curvature, an azimuth at the perigee
  /// or a @c z0: all three of those come from treating the beamline as a third
  /// point on the circle, which is exactly the assumption large radius
  /// tracking drops. So the doublet stage keeps to the non-bending plane and
  /// everything else waits for a triplet fit, which is also where a strip end
  /// is resolved along its strip, the fitted tangent being the first direction
  /// good enough to do it with.
  ///
  /// An edge therefore carries the triplets it is the inner edge of rather
  /// than parameters of its own, and two edges are matched through those. The
  /// levels and chain heads it produces are the prompt graph's, so the rest of
  /// the workflow is unchanged.
  class DisplacedGbtsGraph{

    public:

    /// The edge this graph builds, so that the seeder and the filter can
    /// follow a graph to its edge type rather than being told both.
    using EdgeType = detail::DisplacedGbtsEdge;

    /// Configuration, shared with the prompt graph so that a caller can hand
    /// the same object to either.
    using Config = GbtsGraphConfig;

    /// Takes the same objects as the prompt graph, differing only in the
    /// configuration it is given.
    /// @param config Configuration for the graph
    /// @param geometry GBTS geometry
    /// @param logger Logging instance
    DisplacedGbtsGraph(const Config& config, std::shared_ptr<const GbtsGeometry> geometry,
                       std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
                           "DisplacedGbtsGraph", Acts::Logging::Level::INFO));

    /// Access the configuration, which also carries the chain selection that
    /// seed extraction has to agree with.
    /// @return The configuration
    const Config& config() const { return m_cfg; };

    /// Build the displaced doublet graph from nodes.
    /// @param roi Region of interest descriptor
    /// @param nodeStorage Data storage containing nodes
    /// @param edgeStorage Storage for generated edges
    /// @param bFieldInZ Magnetic field in z, in GeV/(e*mm)
    /// @return Pair of edge count and edge link count
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

    /// Fit the circle through the three nodes of a triplet.
    ///
    /// The transverse plane is inverted about the middle node, which maps
    /// every circle through that node onto a straight line, so the fit is the
    /// line through the two remaining nodes and needs no iteration.
    ///
    /// @param points The three node positions, inside out
    /// @return The circle, or nothing if the nodes are degenerate
    std::optional<detail::TripletCircle> fitTripletCircle(
        const std::array<std::array<float, 3>, 3>& points) const;

    Config m_cfg;

    std::shared_ptr<const GbtsGeometry> m_geometry;

    std::unique_ptr<const Acts::Logger> m_logger;

    const Acts::Logger& logger() const { return *m_logger; }
  };

  
}