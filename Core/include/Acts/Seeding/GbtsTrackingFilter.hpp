// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style
#include "Acts/Seeding/GbtsDataStorage.hpp"  //includes geo which has trigindetsilayer, may move this to trigbase
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <list>
#include <vector>

namespace Acts::Experimental {

struct TrigFTF_GNN_EdgeState {

public:

struct Compare {
    bool operator()(const struct TrigFTF_GNN_EdgeState* s1, const struct TrigFTF_GNN_EdgeState* s2) {
      return s1->m_J > s2->m_J;
    }
  };


  TrigFTF_GNN_EdgeState() {};

  TrigFTF_GNN_EdgeState(bool f) : m_initialized(f) {};

  ~TrigFTF_GNN_EdgeState() {};

  void initialize(TrigFTF_GNN_Edge*);
  void clone(const struct TrigFTF_GNN_EdgeState&);

  float m_J{};

  std::vector<TrigFTF_GNN_Edge*> m_vs;

  float m_X[3]{}, m_Y[2]{}, m_Cx[3][3]{}, m_Cy[2][2]{};
  float m_refX{}, m_refY{}, m_c{}, m_s{};
  
  bool m_initialized{false};

};

#define MAX_EDGE_STATE 2500

class TrigFTF_GNN_TrackingFilter {
 public:
  TrigFTF_GNN_TrackingFilter(const std::vector<TrigInDetSiLayer>&, std::vector<TrigFTF_GNN_Edge>&);
  ~TrigFTF_GNN_TrackingFilter(){};

  void followTrack(TrigFTF_GNN_Edge*, TrigFTF_GNN_EdgeState&);

 protected:

  void propagate(TrigFTF_GNN_Edge*, TrigFTF_GNN_EdgeState&);

  bool update(TrigFTF_GNN_Edge*, TrigFTF_GNN_EdgeState&);

  int getLayerType(int);  


  const std::vector<TrigInDetSiLayer>& m_geo;
  
  std::vector<TrigFTF_GNN_Edge>& m_segStore;
 
  std::vector<TrigFTF_GNN_EdgeState*> m_stateVec;

  TrigFTF_GNN_EdgeState m_stateStore[MAX_EDGE_STATE];

  int m_globalStateCounter{0};

};



}  // namespace Acts::Experimental
