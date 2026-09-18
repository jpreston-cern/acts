#pragma once 

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