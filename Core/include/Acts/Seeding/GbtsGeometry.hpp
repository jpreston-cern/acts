// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style
#include "Acts/TrackFinding/GbtsConnector.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

namespace Acts::Experimental {
class TrigInDetSiLayer {
 public:
  int m_subdet;  // combined ID
  int m_type;    // 0: barrel, +/-n : endcap
  float m_refCoord;
  float m_minBound, m_maxBound;

  TrigInDetSiLayer(int subdet, short int type, float center, float min,
                   float max)
      : m_subdet(subdet),
        m_type(type),
        m_refCoord(center),
        m_minBound(min),
        m_maxBound(max) {}
};

class TrigFTF_GNN_Layer {
public:
  TrigFTF_GNN_Layer(const TrigInDetSiLayer&, float, int);
  ~TrigFTF_GNN_Layer();

  int getEtaBin(float, float) const;

  float getMinBinRadius(int) const;
  float getMaxBinRadius(int) const;

  int num_bins() const {return m_bins.size();} 

  bool verifyBin(const TrigFTF_GNN_Layer*, int, int, float, float) const;

  const TrigInDetSiLayer& m_layer;
  std::vector<int> m_bins;//eta-bin indices
  std::vector<float> m_minRadius;
  std::vector<float> m_maxRadius;
  std::vector<float> m_minBinCoord;
  std::vector<float> m_maxBinCoord;

  float m_minEta, m_maxEta, m_etaBin;

protected:

  float m_etaBinWidth;

  float m_r1, m_z1, m_r2, m_z2;
  int m_nBins; 

};

class TrigFTF_GNN_Geometry {
public:
  TrigFTF_GNN_Geometry(const std::vector<TrigInDetSiLayer>&, const std::unique_ptr<GNN_FasTrackConnector>&);
  ~TrigFTF_GNN_Geometry();
  
  const TrigFTF_GNN_Layer* getTrigFTF_GNN_LayerByKey(unsigned int) const;
  const TrigFTF_GNN_Layer* getTrigFTF_GNN_LayerByIndex(int) const;

  int num_bins() const {return m_nEtaBins;}
  unsigned int num_layers() const {return m_layArray.size();}
  const std::vector<std::pair<int, std::vector<int> > >& bin_groups() const {return m_binGroups;}
  
protected:

  const TrigFTF_GNN_Layer* addNewLayer(const TrigInDetSiLayer&, int);

  float m_etaBinWidth;

  std::map<unsigned int, TrigFTF_GNN_Layer*> m_layMap;
  std::vector<TrigFTF_GNN_Layer*> m_layArray;
  
  int m_nEtaBins;

  std::vector<std::pair<int, std::vector<int> > > m_binGroups;
  
};

}  // namespace Acts::Experimental
