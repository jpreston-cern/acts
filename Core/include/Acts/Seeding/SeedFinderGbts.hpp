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

class SeedingToolBase{
 public:

  SeedingToolBase(const SeedFinderGbtsConfig& config,
                  const TrigFTF_GNN_Geometry* gbtsGeo,
	                const std::vector<TrigInDetSiLayer>& layerGeometry,
                  std::unique_ptr<const Acts::Logger> logger);

  typedef TrigFTF_GNN_Node GNN_Node;
  typedef TrigFTF_GNN_DataStorage GNN_DataStorage;
  typedef TrigFTF_GNN_Edge GNN_Edge;

  SeedContainer2 CreateSeeds(
	const RoiDescriptor& roi, const auto& SpContainerComponents, 
	int max_layers);
  
  
  std::vector<std::vector<SeedingToolBase::GNN_Node>> 
    CreateNodes(const auto& container, int MaxLayers);

  std::pair<int, int> buildTheGraph(const RoiDescriptor&, const std::unique_ptr<GNN_DataStorage>&, std::vector<GNN_Edge>&) const;

  int runCCA(int, std::vector<GNN_Edge>&) const;

 private: 
  
  SeedFinderGbtsConfig m_config;

  std::unique_ptr<GNN_DataStorage> m_storage = nullptr;

  std::vector<TrigInDetSiLayer> m_layerGeometry;

  std::unique_ptr<const Acts::Logger> m_logger =
      Acts::getDefaultLogger("Finder", Acts::Logging::Level::INFO);

  const Acts::Logger &logger() const { return *m_logger; }

};

}  // namespace Acts::Experimental

#include "Acts/Seeding/SeedFinderGbts.ipp"
