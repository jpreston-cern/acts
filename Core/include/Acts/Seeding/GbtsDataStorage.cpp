/*
  Copyright (C) 2002-2025 CERN for the benefit of the ATLAS collaboration
*/

#include "GNN_DataStorage.h"

#include<cmath>
#include<cstring>
#include<algorithm>
namespace Acts::Experimental {

TrigFTF_GNN_EtaBin::TrigFTF_GNN_EtaBin(): m_minRadius(0), m_maxRadius(0) {

  m_in.clear();
  m_vn.clear();
  m_params.clear();
  m_vn.reserve(1000);
}

TrigFTF_GNN_EtaBin::~TrigFTF_GNN_EtaBin() {
  m_in.clear();
  m_vn.clear();
  m_params.clear();
}

void TrigFTF_GNN_EtaBin::sortByPhi() {
  
  std::vector<std::pair<float, const TrigFTF_GNN_Node*> > phiBuckets[32];

  int nBuckets = 31;

  for(const auto& n : m_vn) {
    int bIdx = (int)(0.5*nBuckets*(n->phi()/(float)M_PI + 1.0f));
    phiBuckets[bIdx].push_back(std::make_pair(n->phi(), n));
  }

  for(auto& b : phiBuckets) {
    std::sort(b.begin(), b.end());
  }

  int idx = 0;
  for(const auto& b : phiBuckets) {
    for(const auto& p : b) {
      m_vn[idx++] = p.second;
    }
  }

}


void TrigFTF_GNN_EtaBin::initializeNodes() {
  
  if(m_vn.empty()) return;
  
  m_params.resize(m_vn.size());
  
  m_in.resize(m_vn.size());
  for(auto& v : m_in) v.reserve(50);//reasonably high number of incoming edges per node
  
  std::transform(m_vn.begin(), m_vn.end(), m_params.begin(),
                   [](const TrigFTF_GNN_Node* pN) { std::array<float,5> a = {-100.0, 100.0, pN->phi(), pN->r(), pN->z()}; return a;});
    
  auto [min_iter, max_iter] = std::minmax_element(m_vn.begin(), m_vn.end(),
						  [](const TrigFTF_GNN_Node* s, const TrigFTF_GNN_Node* s1) { return (s->r() < s1->r()); });
  m_maxRadius = (*max_iter)->r();
  m_minRadius = (*min_iter)->r();
}

void TrigFTF_GNN_EtaBin::generatePhiIndexing(float dphi) {

  for(unsigned int nIdx=0;nIdx<m_vn.size();nIdx++) {

    float phi = m_params[nIdx][2];
    if(phi <= M_PI-dphi) continue;    
    m_vPhiNodes.push_back(std::pair<float, unsigned int>(phi - 2*M_PI, nIdx));
    
  }

  for(unsigned int nIdx=0;nIdx<m_vn.size();nIdx++) {
    float phi = m_params[nIdx][2];
    m_vPhiNodes.push_back(std::pair<float, unsigned int>(phi, nIdx));
  }

  for(unsigned int nIdx=0;nIdx<m_vn.size();nIdx++) {

    float phi = m_params[nIdx][2];
    if(phi >= -M_PI + dphi) break;
    m_vPhiNodes.push_back(std::pair<float, unsigned int>(phi + 2*M_PI, nIdx));
  }
  
}

TrigFTF_GNN_DataStorage::TrigFTF_GNN_DataStorage(const TrigFTF_GNN_Geometry& g) : m_geo(g) {
  m_etaBins.resize(g.num_bins());
}


TrigFTF_GNN_DataStorage::~TrigFTF_GNN_DataStorage() {

}

int TrigFTF_GNN_DataStorage::loadPixelGraphNodes(short layerIndex, const std::vector<TrigFTF_GNN_Node>& coll, bool useML) {

  int nLoaded = 0;

  const TrigFTF_GNN_Layer* pL = m_geo.getTrigFTF_GNN_LayerByIndex(layerIndex);

  if(pL == nullptr) {
    return -1;
  }
  
  bool isBarrel = (pL->m_layer.m_type == 0);
 
  for(const auto& node : coll) {

    int binIndex = pL->getEtaBin(node.z(), node.r());

    if(binIndex == -1) {
      continue;
    }
    
    if(isBarrel) {
      m_etaBins.at(binIndex).m_vn.push_back(&node);
    }
    else {
      if (useML) {
	float cluster_width = node.pixelClusterWidth();
	if(cluster_width > 0.2) continue;
      }
      m_etaBins.at(binIndex).m_vn.push_back(&node);
    }
    
    nLoaded++;
    
  }
  
  return nLoaded;
}

int TrigFTF_GNN_DataStorage::loadStripGraphNodes(short layerIndex, const std::vector<TrigFTF_GNN_Node>& coll) {

  int nLoaded = 0;

  const TrigFTF_GNN_Layer* pL = m_geo.getTrigFTF_GNN_LayerByIndex(layerIndex);

  if(pL == nullptr) {
    return -1;
  }
   
  for(const auto& node : coll) {

    int binIndex = pL->getEtaBin(node.z(), node.r());

    if(binIndex == -1) {
      continue;
    }

    m_etaBins.at(binIndex).m_vn.push_back(&node);
    nLoaded++;
  }
  
  return nLoaded;
}

unsigned int TrigFTF_GNN_DataStorage::numberOfNodes() const {

  unsigned int n=0;
  
  for(const auto& b : m_etaBins) {
    n += b.m_vn.size();
  }
  return n;
}

void TrigFTF_GNN_DataStorage::sortByPhi() {
    
  for(auto& b : m_etaBins) b.sortByPhi();
}

void TrigFTF_GNN_DataStorage::initializeNodes(bool useML) {
  
  for(auto& b : m_etaBins) {
    b.initializeNodes();
  }
  
  if(!useML) return;
  
  unsigned int nL = m_geo.num_layers();

  for(unsigned int layerIdx=0;layerIdx<nL;layerIdx++) {

    const TrigFTF_GNN_Layer* pL = m_geo.getTrigFTF_GNN_LayerByIndex(layerIdx);

    if(pL->m_layer.m_subdet < 20000) {//skip strips volumes: layers in range [1200X-1400X]
      continue;
    }
    
    bool isBarrel = (pL->m_layer.m_type == 0);

    if(!isBarrel) continue;

    int nBins = pL->m_bins.size();

    for(int b=0;b<nBins;b++) {//loop over eta-bins in Layer

      TrigFTF_GNN_EtaBin& B = m_etaBins.at(pL->m_bins.at(b));

      if(B.empty()) continue;
      
      for(unsigned int nIdx=0;nIdx<B.m_vn.size();nIdx++) {
        float cluster_width = B.m_vn[nIdx]->pixelClusterWidth();
	//adjusting cuts using fitted boundaries of |cot(theta)| vs. cluster z-width distribution 
        float min_tau = 6.7*(cluster_width - 0.2);//linear fit
        float max_tau = 1.6 + 0.15/(cluster_width + 0.2) + 6.1*(cluster_width - 0.2);//linear fit + correction for short clusters

        B.m_params[nIdx][0] = min_tau;
        B.m_params[nIdx][1] = max_tau;
        
      }
    }
  }
}

void TrigFTF_GNN_DataStorage::generatePhiIndexing(float dphi) {
  for(auto& b : m_etaBins) b.generatePhiIndexing(dphi);
}

}
