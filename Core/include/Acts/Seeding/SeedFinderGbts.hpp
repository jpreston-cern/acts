// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Seed.hpp"
#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"
#include "Acts/TrackFinding/RoiDescriptor.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "Acts/TrackFinding/GbtsConnector.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsDataStorage.hpp"

#include <memory>
#include <string>
#include <vector>

namespace Acts::Experimental {
  
//defined the tuple template used to carry the spacepoint components
using SPContainerComponentsType = std::tuple<
                                             SpacePointContainer2,
                                             SpacePointColumnProxy<int, true>,
                                             SpacePointColumnProxy<float, true>>;

class SeedFinderGbts{
 public:

  SeedFinderGbts( const SeedFinderGbtsConfig config,
                   const GbtsGeometry* gbtsGeo,
	                 const std::vector<TrigInDetSiLayer>* layerGeometry,
                   std::unique_ptr<const Acts::Logger> logger);

  typedef GbtsNode GNN_Node;
  typedef GbtsDataStorage GNN_DataStorage;
  typedef GbtsEdge GNN_Edge;

  SeedContainer2 CreateSeeds(
	  const RoiDescriptor& roi, 
    const SPContainerComponentsType& SpContainerComponents, 
	  int max_layers);
  
  
  std::vector<std::vector<SeedFinderGbts::GNN_Node>> 
    CreateNodes(const auto& container, int MaxLayers);

  std::pair<int, int> buildTheGraph(const RoiDescriptor&, const std::unique_ptr<GNN_DataStorage>&, std::vector<GNN_Edge>&) const;

  int runCCA(int, std::vector<GNN_Edge>&) const;

 private: 
  
  SeedFinderGbtsConfig m_config;

  const GbtsGeometry* m_geo;

  std::unique_ptr<GNN_DataStorage> m_storage = nullptr;

  const std::vector<TrigInDetSiLayer>* m_layerGeometry;

  std::unique_ptr<const Acts::Logger> m_logger =
      Acts::getDefaultLogger("Finder", Acts::Logging::Level::INFO);

  const Acts::Logger &logger() const { return *m_logger; }

};

}  // namespace Acts::Experimental


