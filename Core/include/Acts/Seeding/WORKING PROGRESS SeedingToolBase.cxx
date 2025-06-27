/*
  Copyright (C) 2002-2025 CERN for the benefit of the ATLAS collaboration
*/

#include "InDetIdentifier/SCT_ID.h"
#include "InDetIdentifier/PixelID.h" 

#include "AtlasDetDescr/AtlasDetectorID.h"

#include "PathResolver/PathResolver.h"

#include "GNN_TrackingFilter.h"

#include "IRegionSelector/IRegSelTool.h"

#include "SeedingToolBase.h"


std::pair<int, int> SeedingToolBase::buildTheGraph(const IRoiDescriptor& roi, const std::unique_ptr<TrigFTF_GNN_DataStorage>& storage, std::vector<TrigFTF_GNN_Edge>& edgeStorage) const {

  const float M_2PI = 2.0*M_PI;
  
  const float cut_dphi_max      = m_LRTmode ? 0.07 : 0.012;
  const float cut_dcurv_max     = m_LRTmode ? 0.015 : 0.001;
  const float cut_tau_ratio_max = m_LRTmode ? 0.015 : 0.007;
  const float min_z0            = m_LRTmode ? -600.0 : roi.zedMinus();
  const float max_z0            = m_LRTmode ? 600.0 : roi.zedPlus();
  const float min_deltaPhi      = m_LRTmode ? 0.01f : 0.001f;
  
  const float maxOuterRadius    = m_LRTmode ? 1050.0 : 550.0;

  const float cut_zMinU = min_z0 + maxOuterRadius*roi.dzdrMinus();
  const float cut_zMaxU = max_z0 + maxOuterRadius*roi.dzdrPlus();

  const float ptCoeff = 0.29997*1.9972/2.0;// ~0.3*B/2 - assuming nominal field of 2*T

  float tripletPtMin = 0.8*m_minPt;//correction due to limited pT resolution
  
  float maxCurv = ptCoeff/tripletPtMin;
 
  const float maxKappa_high_eta          = m_LRTmode ? 1.0*maxCurv : std::sqrt(0.8)*maxCurv;
  const float maxKappa_low_eta           = m_LRTmode ? 1.0*maxCurv : std::sqrt(0.6)*maxCurv;
  const float dphi_coeff                 = m_LRTmode ? 1.0*maxCurv : 0.68*maxCurv;
  
  const float minDeltaRadius = 2.0;
    
  float deltaPhi = 0.5f*m_phiSliceWidth;//the default sliding window along phi
 
  unsigned int nConnections = 0;
  
  edgeStorage.reserve(m_nMaxEdges);
  
  int nEdges = 0;

  for(const auto& bg : m_geo->bin_groups()) {//loop over bin groups
    
    TrigFTF_GNN_EtaBin& B1 = storage->getEtaBin(bg.first);

    if(B1.empty()) continue;

    float rb1 = B1.getMinBinRadius();
 
    for(const auto& b2_idx : bg.second) {

      const TrigFTF_GNN_EtaBin& B2 = storage->getEtaBin(b2_idx);

      if(B2.empty()) continue;
      
      float rb2 = B2.getMaxBinRadius();
    
      if(m_useEtaBinning) {
	deltaPhi = min_deltaPhi + dphi_coeff*std::fabs(rb2-rb1);	
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
		
	  if (m_doubletFilterRZ) {
		  
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
	  
	  if (m_matchBeforeCreate) {//match edge candidate against edges incoming to n2

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
	
	  if(nEdges < m_nMaxEdges) {
	  
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

