// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/SeedFinderGbts.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp" 

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

namespace Acts::Experimental {


SeedingToolBase::SeedingToolBase(
	SeedFinderGbtsConfig config,
    const TrigFTF_GNN_Geometry* gbtsGeo,
	const std::vector<TrigInDetSiLayer>* layerGeometry,
    std::unique_ptr<const Acts::Logger> logger) 
    : m_config(config),
	  m_geo(gbtsGeo),
      m_storage(std::make_unique<GNN_DataStorage>(*m_geo)),
	  m_layerGeometry(layerGeometry),
	  m_logger(std::move(logger)) {}


SeedContainer2 SeedingToolBase::CreateSeeds(
	const RoiDescriptor& roi, const auto& SpContainerComponents, 
	int max_layers){

	SeedContainer2 SeedContainer;

	std::vector<std::vector<GNN_Node>> node_storage =
	 CreateNodes(SpContainerComponents, max_layers);
	
	unsigned int nPixelLoaded = 0;
    unsigned int nStripLoaded = 0;

	//load the nodes into storage
	for(size_t l = 0; l < node_storage.size(); l++) {

      const std::vector<GNN_Node>& nodes = node_storage[l];

      if(nodes.size() == 0) continue;

	  bool is_pixel = true;
      if(is_pixel) //placeholder for now until strip hits are added in
		nPixelLoaded += m_storage->loadPixelGraphNodes(l, nodes, m_config.m_useML);
      else
		nStripLoaded += m_storage->loadStripGraphNodes(l, nodes);
    }

	ACTS_DEBUG("Loaded "<<nPixelLoaded<<" pixel spacepoints and "<<nStripLoaded<<" strip spacepoints");

    m_storage->sortByPhi();

    m_storage->initializeNodes(m_config.m_useML);

	m_config.m_phiSliceWidth = 2*std::numbers::pi/m_config.m_nMaxPhiSlice;
    m_storage->generatePhiIndexing(1.5*m_config.m_phiSliceWidth);

	std::vector<GNN_Edge> edgeStorage; 

    std::pair<int, int> graphStats = buildTheGraph(roi, m_storage, edgeStorage);

    ACTS_DEBUG("Created graph with "<<graphStats.first<<" edges and "<<graphStats.second<< " edge links");

	int maxLevel = runCCA(graphStats.first, edgeStorage);

    ACTS_DEBUG("Reached Level "<<maxLevel<<" after GNN iterations");

    int minLevel = 3;//a triplet + 2 confirmation

    if(m_config.m_LRTmode) {
        minLevel = 2;//a triplet + 1 confirmation
    }

    

	std::vector<GNN_Edge*> vSeeds;

    vSeeds.reserve(graphStats.first/2);

    for(int edgeIndex = 0; edgeIndex < graphStats.first; edgeIndex++) {
      
      GNN_Edge* pS = &(edgeStorage.at(edgeIndex));

      if(pS->m_level < minLevel) continue;

      vSeeds.push_back(pS);
    }

 
    std::sort(vSeeds.begin(), vSeeds.end(), GNN_Edge::CompareLevel());

    //backtracking

    TrigFTF_GNN_TrackingFilter tFilter(*m_layerGeometry, edgeStorage); //add this file in
	
	for(auto pS : vSeeds) {
        if(pS->m_level == -1) continue;

        TrigFTF_GNN_EdgeState rs(false); //this is in tracking filter as well

        tFilter.followTrack(pS, rs);

        if(!rs.m_initialized) {
            continue;
        }

        if(static_cast<int>(rs.m_vs.size()) < minLevel) continue;

        std::vector<const GNN_Node*> vN;

        for(std::vector<GNN_Edge*>::reverse_iterator sIt=rs.m_vs.rbegin();sIt!=rs.m_vs.rend();++sIt) {
            (*sIt)->m_level = -1;//mark as collected
            
            if(sIt == rs.m_vs.rbegin()) {
                vN.push_back((*sIt)->m_n1);
            }

            vN.push_back((*sIt)->m_n2);
	    
        }

        if(vN.size()<3) continue; 
	
       
		//add to seed container:
		std::array<SpacePointIndex2, 3> Sp_Indexes{};
		int index = 0;
		for(const auto* vNptr : vN){

			Sp_Indexes[index] = vNptr->sp_idx();
			index ++;
		}
		
		SeedContainer.createSeed(Sp_Indexes);
    }

    ACTS_DEBUG("GBTS created "<<SeedContainer.size()<<" seeds");

    return SeedContainer;
		
}

std::vector<std::vector<SeedingToolBase::GNN_Node>> 
  SeedingToolBase::CreateNodes(const auto& container, int MaxLayers){

	std::vector<std::vector<SeedingToolBase::GNN_Node>> node_storage;
	//reserve for better efficiency 
	node_storage.resize(MaxLayers);
	for(auto& v : node_storage) v.reserve(100000);

		for(auto sp : std::get<0>(container)){

			//for every sp in container,
			//add its variables to node_storage organised by layer 
			int layer = sp.extra(std::get<1>(container));

			//add node to storage 
			SeedingToolBase::GNN_Node& node = node_storage[layer].emplace_back(layer);

			//fill the node with spacepoint variables
			node.m_x = sp.x();
			node.m_y = sp.y();
			node.m_z = sp.z();
			node.m_r = sp.r();
			node.m_phi = sp.phi();
			node.m_idx = sp.index();//change node so that is uses SpacePointIndex2 (doesnt affect code but i guess it looks cleaner)
			node.m_pcw = sp.extra(std::get<2>(container));
		}

	return node_storage;

}

std::pair<int, int> SeedingToolBase::buildTheGraph(const RoiDescriptor& roi, const std::unique_ptr<TrigFTF_GNN_DataStorage>& storage, std::vector<TrigFTF_GNN_Edge>& edgeStorage) const {

  const float M_2PI = 2.0*M_PI;
  
  const float cut_dphi_max      = m_config.m_LRTmode ? 0.07f : 0.012f;
  const float cut_dcurv_max     = m_config.m_LRTmode ? 0.015f : 0.001f;
  const float cut_tau_ratio_max = m_config.m_LRTmode ? 0.015f : static_cast<float>(m_config.m_tau_ratio_cut);
  const float min_z0            = m_config.m_LRTmode ? -600.0 : roi.zedMinus();
  const float max_z0            = m_config.m_LRTmode ? 600.0 : roi.zedPlus();
  const float min_deltaPhi      = m_config.m_LRTmode ? 0.01f : 0.001f;
  
  const float maxOuterRadius    = m_config.m_LRTmode ? 1050.0 : 550.0;

  const float cut_zMinU = min_z0 + maxOuterRadius*roi.dzdrMinus();
  const float cut_zMaxU = max_z0 + maxOuterRadius*roi.dzdrPlus();

  const float ptCoeff = 0.29997*1.9972/2.0;// ~0.3*B/2 - assuming nominal field of 2*T

  float tripletPtMin = 0.8*m_config.m_minPt;//correction due to limited pT resolution
  const float pt_scale     = 900.0/m_config.m_minPt;//to re-scale original tunings done for the 900 MeV pT cut

  float maxCurv = ptCoeff/tripletPtMin;
 
  float maxKappa_high_eta          = m_config.m_LRTmode ? 1.0*maxCurv : std::sqrt(0.8)*maxCurv;
  float maxKappa_low_eta           = m_config.m_LRTmode ? 1.0*maxCurv : std::sqrt(0.6)*maxCurv;

  if(!m_config.m_useOldTunings && !m_config.m_LRTmode) {//new settings for curvature cuts
    maxKappa_high_eta          = 4.75e-4f*pt_scale;
    maxKappa_low_eta           = 3.75e-4f*pt_scale;
  }

  const float dphi_coeff                 = m_config.m_LRTmode ? 1.0*maxCurv : 0.68*maxCurv;
  
  const float minDeltaRadius = 2.0;
    
  float deltaPhi = 0.5f*m_config.m_phiSliceWidth;//the default sliding window along phi
 
  unsigned int nConnections = 0;
  
  edgeStorage.reserve(m_config.m_nMaxEdges);
  
  int nEdges = 0;

  for(const auto& bg : m_geo->bin_groups()) {//loop over bin groups
    
    TrigFTF_GNN_EtaBin& B1 = storage->getEtaBin(bg.first);

    if(B1.empty()) continue;

    float rb1 = B1.getMinBinRadius();
 
    for(const auto& b2_idx : bg.second) {

      const TrigFTF_GNN_EtaBin& B2 = storage->getEtaBin(b2_idx);

      if(B2.empty()) continue;
      
      float rb2 = B2.getMaxBinRadius();
    
      if(m_config.m_useEtaBinning) {
		float abs_dr = std::fabs(rb2-rb1);
		if (m_config.m_useOldTunings) {
	  		deltaPhi = min_deltaPhi + dphi_coeff*abs_dr;
		}
		else {
	  		if(abs_dr < 60.0) {
	    deltaPhi = 0.002f + 4.33e-4f*pt_scale*abs_dr;
	  	} 
	  	else {
	    	deltaPhi = 0.015f + 2.2e-4f*pt_scale*abs_dr;
	  	}
	}

	}

      unsigned int first_it = 0;

      for(unsigned int n1Idx = 0;n1Idx<B1.m_vn.size();n1Idx++) {//loop over nodes in Layer 1

	std::vector<unsigned int>& v1In = B1.m_in[n1Idx];   

	if(v1In.size() >= MAX_SEG_PER_NODE) continue;
      
	const std::array<float, 5>& n1pars = B1.m_params[n1Idx];

	float phi1 = n1pars[2];
	float r1 = n1pars[3];
	float z1 = n1pars[4];
      
	//sliding window phi1 +/- deltaPhi
      
	float minPhi = phi1 - deltaPhi;
	float maxPhi = phi1 + deltaPhi;
      
	for(unsigned int n2PhiIdx = first_it; n2PhiIdx<B2.m_vPhiNodes.size();n2PhiIdx++) {//sliding window over nodes in Layer 2
	
	  float phi2 = B2.m_vPhiNodes[n2PhiIdx].first;
	
	  if(phi2 < minPhi) {
	    first_it = n2PhiIdx;
	    continue;
	  }
	  if(phi2 > maxPhi) break;
	
	  unsigned int n2Idx = B2.m_vPhiNodes[n2PhiIdx].second;
	
	  const std::vector<unsigned int>& v2In = B2.m_in[n2Idx];
        
	  if(v2In.size() >= MAX_SEG_PER_NODE) continue;
		
	  const std::array<float, 5>& n2pars = B2.m_params[n2Idx];
	
	  float r2 = n2pars[3];
	  
	  float dr = r2 - r1;
	
	  if(dr < minDeltaRadius) {
	    continue;
	  }
	
	  float z2 = n2pars[4];

	  float dz = z2 - z1;
	  float tau = dz/dr;
	  float ftau = std::fabs(tau);
	  if (ftau > 36.0) {
	    continue;
	  }
	
	  if(ftau < n1pars[0]) continue;
	  if(ftau > n1pars[1]) continue;

	  if(ftau < n2pars[0]) continue;
	  if(ftau > n2pars[1]) continue;
		
	  if (m_config.m_doubletFilterRZ) {
		  
	    float z0 = z1 - r1*tau;
	  
	    if(z0 < min_z0 || z0 > max_z0) continue;
	  
	    float zouter = z0 + maxOuterRadius*tau;
	  
	    if(zouter < cut_zMinU || zouter > cut_zMaxU) continue;                
	  }
		
	  float curv = (phi2-phi1)/dr;
	  float abs_curv = std::abs(curv);
		
	  if(ftau < 4.0) {//eta = 2.1
	    if(abs_curv > maxKappa_low_eta) {
	      continue;
	    }
	  }
	  else {
	    if(abs_curv > maxKappa_high_eta) {
	      continue;
	    }
	  }

	  float exp_eta = std::sqrt(1+tau*tau)-tau;
	  
	  if (m_config.m_matchBeforeCreate) {//match edge candidate against edges incoming to n2

	    bool isGood = v2In.size() <= 2;//we must have enough incoming edges to decide

	    if(!isGood) {

	      float uat_1 = 1.0f/exp_eta;
		    
	      for(const auto& n2_in_idx : v2In) {
		    
		float tau2 = edgeStorage.at(n2_in_idx).m_p[0]; 
		float tau_ratio = tau2*uat_1 - 1.0f;
		
		if(std::fabs(tau_ratio) > cut_tau_ratio_max){//bad match
		  continue;
		}
		isGood = true;//good match found
		break;
	      }
	    }
	    
	    if(!isGood) {//no match found, skip creating [n1 <- n2] edge
	      continue;
	    }
	  }
	  
	  float dPhi2 = curv*r2;
	  float dPhi1 = curv*r1;
	
	  if(nEdges < m_config.m_nMaxEdges) {
	  
	    edgeStorage.emplace_back(B1.m_vn[n1Idx], B2.m_vn[n2Idx], exp_eta, curv, phi1 + dPhi1);
	    
	    if(v1In.size() < MAX_SEG_PER_NODE) v1In.push_back(nEdges);
		  
	    int outEdgeIdx = nEdges;
	  
	    float uat_2  = 1/exp_eta;
	    float Phi2  = phi2 + dPhi2;
	    float curv2 = curv;
	    
	    for(const auto& inEdgeIdx : v2In) {//looking for neighbours of the new edge
	    
	      TrigFTF_GNN_Edge* pS = &(edgeStorage.at(inEdgeIdx));
	      
	      if(pS->m_nNei >= N_SEG_CONNS) continue;
	    
	      float tau_ratio = pS->m_p[0]*uat_2 - 1.0f;
	      
	      if(std::abs(tau_ratio) > cut_tau_ratio_max){//bad match
		continue;
	      }
	      
	      float dPhi =  Phi2 - pS->m_p[2];
	      
	      if(dPhi<-M_PI) dPhi += M_2PI;
	      else if(dPhi>M_PI) dPhi -= M_2PI;
	      
	      if(dPhi < -cut_dphi_max || dPhi > cut_dphi_max) {
		continue;
	      }
            
	      float dcurv = curv2 - pS->m_p[1];
            
	      if(dcurv < -cut_dcurv_max || dcurv > cut_dcurv_max) {
		continue;
	      }
            
	      pS->m_vNei[pS->m_nNei++] = outEdgeIdx;
	    
	      nConnections++;
	    
	    }
	    nEdges++;		
	  }
	} //loop over n2 (outer) nodes
      } //loop over n1 (inner) nodes
    } //loop over bins in Layer 2
  } //loop over bin groups

  if(nEdges >= m_config.m_nMaxEdges) {
    ACTS_WARNING("Maximum number of graph edges exceeded - possible efficiency loss "<< nEdges);
  }
  return std::make_pair(nEdges, nConnections);
}

int SeedingToolBase::runCCA(int nEdges, std::vector<TrigFTF_GNN_Edge>& edgeStorage) const {

  const int maxIter = 15;

  int maxLevel = 0;

  int iter = 0;
  
  std::vector<TrigFTF_GNN_Edge*> v_old;
  
  for(int edgeIndex=0;edgeIndex<nEdges;edgeIndex++) {

    TrigFTF_GNN_Edge* pS = &(edgeStorage[edgeIndex]);
    if(pS->m_nNei == 0) continue;
    
    v_old.push_back(pS);//TO-DO: increment level for segments as they already have at least one neighbour
  }

  for(;iter<maxIter;iter++) {

    //generate proposals
    std::vector<TrigFTF_GNN_Edge*> v_new;
    v_new.clear();
    v_new.reserve(v_old.size());
    
    for(auto pS : v_old) {
      
      int next_level = pS->m_level;
          
      for(int nIdx=0;nIdx<pS->m_nNei;nIdx++) {
	
        unsigned int nextEdgeIdx = pS->m_vNei[nIdx];
            
        TrigFTF_GNN_Edge* pN = &(edgeStorage[nextEdgeIdx]);
            
        if(pS->m_level == pN->m_level) {
          next_level = pS->m_level + 1;
          v_new.push_back(pS);
          break;
        }
      }
      
      pS->m_next = next_level;//proposal
    }
  
    //update

    int nChanges = 0;
      
    for(auto pS : v_new) {
      if(pS->m_next != pS->m_level) {
        nChanges++;
        pS->m_level = pS->m_next;
        if(maxLevel < pS->m_level) maxLevel = pS->m_level;
      }
    }

    if(nChanges == 0) break;


    v_old = std::move(v_new);
    v_new.clear();
  }

  return maxLevel;  
}



}  // namespace Acts::Experimental
