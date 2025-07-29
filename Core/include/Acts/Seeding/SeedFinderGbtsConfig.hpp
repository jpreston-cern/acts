// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include <memory>

// core algorithm so in acts namespace
namespace Acts::Experimental {

struct SeedFinderGbtsConfig {
  
 // GbtsSeedingAlgorithm options 
  bool BeamSpotCorrection = false;
  std::string ConnectorInputFile{};  // Path to the connector configuration file
                                     // that defines the layer connections

 
 //SeedFinderGbts option
  bool m_LRTmode = false;
  bool m_useML = false; //use cluster width 
  bool m_matchBeforeCreate = false;
  bool m_useOldTunings = false;
  bool m_tau_ratio_cut = 0.007;
  float m_etaBinOverride = 0.0f; //specify non-zero to override eta bin width from connection file (default 0.2 in createLinkingScheme.py)
  float m_nMaxPhiSlice = 53; // used to calculate phi slices
  float m_minPt = 1000.0;
  float m_phiSliceWidth{}; //derived in CreatSeeds function

    //BuildTheGraph() options
  double ptCoeff = 0.29997 * 1.9972 / 2.0; // ~0.3*B/2 - assumes nominal field of 2*T
  float minPt = 400. * Acts::UnitConstants::MeV;
  bool m_useEtaBinning = true;  // bool to use eta binning from geometry structure
  bool m_doubletFilterRZ = true;  // bool applies new Z cuts on doublets
  int m_nMaxEdges = 2000000;     // max number of Gbts edges/doublets
  float m_minDeltaRadius = 2.0;

  
  
   
  
//unsure if i still need these
  float sigmaScattering = 5;  
  float highland = 0;
  float maxScatteringAngle2 = 0;
  bool m_useClusterWidth = false; //defineitly get rid of this one 
  std::vector<TrigInDetSiLayer> m_layerGeometry; //get rid off

  

  ////
  // 2 member functions
  SeedFinderGbtsConfig calculateDerivedQuantities() const {
    // thorw statement if the isInternalUnits member is false, ie if dont call
    // this function
    SeedFinderGbtsConfig config = *this;
    // use a formula to calculate scattering

    return config;
  }

  SeedFinderGbtsConfig toInternalUnits() const {
    // throw statement if the isInternalUnits member is false, ie if dont call
    // this function
    SeedFinderGbtsConfig config = *this;
    // divides inputs by 1mm, all ones input
    // changes member inInInternalUnits to true
    return config;
  };

};  // end of config struct

}  // namespace Acts::Experimental
