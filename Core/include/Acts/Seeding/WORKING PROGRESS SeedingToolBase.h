/*
  Copyright (C) 2002-2025 CERN for the benefit of the ATLAS collaboration
*/

#pragma once 



#include <string>
#include <vector>

#include "GNN_FasTrackConnector.h"
#include "GNN_Geometry.h"
#include "GNN_DataStorage.h"



class SeedingToolBase{
 public:

  SeedingToolBase(SeedFinderGbtsConfig);
  
 


  typedef TrigFTF_GNN_Node GNN_Node;
  typedef TrigFTF_GNN_DataStorage GNN_DataStorage;
  typedef TrigFTF_GNN_Edge GNN_Edge;

  StatusCode createSeeds(
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


