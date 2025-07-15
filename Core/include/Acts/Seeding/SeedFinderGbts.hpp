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

#include "Acts/Trackfinding/GbtsConnector.h"
#include "Acts/Seeding/GbtsGeometry.h"
#include "Acts/Seeding/GBtsDataStorage.h"

#include <memory>
#include <string>
#include <vector>

namespace Acts::Experimental {

class SeedingToolBase{
 public:

  SeedingToolBase(SeedFinderGbtsConfig);
  
  typedef TrigFTF_GNN_Node GNN_Node;
  typedef TrigFTF_GNN_DataStorage GNN_DataStorage;
  typedef TrigFTF_GNN_Edge GNN_Edge;

  SeedContainer2 createSeeds(
    const RoiDescriptor& roi, const Acts::Experimental::SpacePointContainer2& SPcontainer, 
	  int max_layers, Acts::Experimental::SeedContainer2& SeedContainer);
  
  
  std::vector<std::vector<Acts::Experimental::GNN_Node>> 
    CreateNodes(const Acts::Experimental::SpacePointContainer2& SPcontainer, int MaxLayers);

  std::pair<int, int> buildTheGraph(const IRoiDescriptor&, const std::unique_ptr<GNN_DataStorage>&, std::vector<GNN_Edge>&) const;

  int runCCA(int, std::vector<GNN_Edge>&) const;

 private: 
  
  SeedFinderGbtsConfig m_config;

  std::unique_ptr<GbtsDataStorage<external_spacepoint_t>> m_storage{nullptr};

  std::vector<TrigInDetSiLayer> m_layerGeometry;

};

}  // namespace Acts::Experimental

#include "Acts/Seeding/SeedFinderGbts.ipp"
