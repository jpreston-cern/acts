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
  SeedingToolBase() = default;
  
 


  typedef TrigFTF_GNN_Node GNN_Node;
  typedef TrigFTF_GNN_DataStorage GNN_DataStorage;
  typedef TrigFTF_GNN_Edge GNN_Edge;

  
  
  std::vector<std::vector<Acts::Experimental::GNN_Node>> 
    CreateNodes(Acts::Experimental::SpacePointContainer2& container, int MaxLayers);

  std::pair<int, int> buildTheGraph(const IRoiDescriptor&, const std::unique_ptr<GNN_DataStorage>&, std::vector<GNN_Edge>&) const;

  int runCCA(int, std::vector<GNN_Edge>&) const;

 private: 
  
  
};


