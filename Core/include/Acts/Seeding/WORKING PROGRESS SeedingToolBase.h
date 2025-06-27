/*
  Copyright (C) 2002-2025 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGINDETPATTRECOTOOLS_SEEDINGTOOLBASE_H
#define TRIGINDETPATTRECOTOOLS_SEEDINGTOOLBASE_H

#include "GaudiKernel/ToolHandle.h"
//#include "TrigInDetToolInterfaces/ITrigInDetTrackSeedingTool.h"
#include "AthenaBaseComps/AthAlgTool.h"
#include "StoreGate/ReadHandleKey.h"
#include <string>
#include <vector>

#include "IRegionSelector/IRegSelTool.h"
#include "TrigInDetToolInterfaces/ITrigL2LayerNumberTool.h"

#include "GNN_FasTrackConnector.h"
#include "GNN_Geometry.h"
#include "GNN_DataStorage.h"

class AtlasDetectorID;
class SCT_ID;
class PixelID;

class SeedingToolBase{
 public:
  SeedingToolBase() = default;
  
 protected:

  typedef TrigFTF_GNN_Node GNN_Node;
  typedef TrigFTF_GNN_DataStorage GNN_DataStorage;
  typedef TrigFTF_GNN_Edge GNN_Edge;

  
  

  std::pair<int, int> buildTheGraph(const IRoiDescriptor&, const std::unique_ptr<GNN_DataStorage>&, std::vector<GNN_Edge>&) const;

  int runCCA(int, std::vector<GNN_Edge>&) const;
  
  ToolHandle<ITrigL2LayerNumberTool> m_layerNumberTool {this, "layerNumberTool", "TrigL2LayerNumberToolITk"};


  
};

#endif
