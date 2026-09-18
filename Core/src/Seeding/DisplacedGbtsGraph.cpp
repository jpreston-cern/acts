#include "Acts/Seeding/DisplacedGbtsGraph.hpp"

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

  
} // Acts::Experimental namespace