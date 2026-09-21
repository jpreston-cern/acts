
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

  class DisplacedGbtsGraph{

    public:

    /// The edge this graph builds, so that the seeder and the filter can
    /// follow a graph to its edge type rather than being told both.
    using EdgeType = detail::DisplacedGbtsEdge;

    /// Configuration, shared with the prompt graph so that a caller can hand
    /// the same object to either.
    using Config = GbtsGraphConfig;

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