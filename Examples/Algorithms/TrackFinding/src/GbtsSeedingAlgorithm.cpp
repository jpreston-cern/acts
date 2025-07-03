// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/GbtsSeedingAlgorithm.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <random>
#include <sstream>
#include <vector>

template class Acts::Experimental::GbtsLayer<ActsExamples::SimSpacePoint>; //this is done so that the teplates are initlaised onto the spacepoints during compile (saves time during running)
template class Acts::Experimental::GbtsGeometry<ActsExamples::SimSpacePoint>;
template class Acts::Experimental::GbtsNode<ActsExamples::SimSpacePoint>;
template class Acts::Experimental::GbtsEtaBin<ActsExamples::SimSpacePoint>;
template struct Acts::Experimental::GbtsSP<ActsExamples::SimSpacePoint>;
template class Acts::Experimental::GbtsDataStorage<ActsExamples::SimSpacePoint>;
template class Acts::Experimental::GbtsEdge<ActsExamples::SimSpacePoint>;

// constructor:
ActsExamples::GbtsSeedingAlgorithm::GbtsSeedingAlgorithm(
    ActsExamples::GbtsSeedingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("SeedingAlgorithm", lvl), m_cfg(std::move(cfg)) {

      std::cout<<"Jasper: Instance of GbtsSeedingAlgorithm created (contructor called) from Examples/Algorithms/TrakFinding/src/GbtsSeedingAlgorithm.cpp"<<std::endl;
  // fill config struct
  //m_cfg.layerMappingFile = m_cfg.layerMappingFile;

  m_cfg.seedFinderConfig =
      m_cfg.seedFinderConfig.toInternalUnits().calculateDerivedQuantities();

  m_cfg.seedFinderOptions =
      m_cfg.seedFinderOptions.toInternalUnits().calculateDerivedQuantities(
          m_cfg.seedFinderConfig);

  for (const auto &spName : m_cfg.inputSpacePoints) {//for every collection of spacepoints
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");// check to see if the spacepoint colletion vector is empty (this is a vector of collection not a single collection)
    }

    auto &handle = m_inputSpacePoints.emplace_back(//add a data handle for the spacepointcontainer to the vector of data handles 
        std::make_unique<ReadDataHandle<SimSpacePointContainer>>(//
            this,
            "InputSpacePoints#" + std::to_string(m_inputSpacePoints.size()))); //define the name of the collection
    handle->initialize(spName); //associate the hadle with the collection string 
  }

  m_outputSeeds.initialize(m_cfg.outputSeeds); //associate the data handle with the string name for output seeds

  m_inputClusters.initialize(m_cfg.inputClusters); //associate the data handle with the string name for input clusters

  // map
  m_cfg.ActsGbtsMap = makeActsGbtsMap();
  // input trig vector
  m_cfg.seedFinderConfig.m_layerGeometry = LayerNumbering();

  std::ifstream input_ifstream(
      m_cfg.seedFinderConfig.ConnectorInputFile.c_str(), std::ifstream::in);

  // connector
  std::unique_ptr<Acts::Experimental::GbtsConnector> inputConnector =
      std::make_unique<Acts::Experimental::GbtsConnector>(input_ifstream);

  m_gbtsGeo = std::make_unique<Acts::Experimental::GbtsGeometry<SimSpacePoint>>(
      m_cfg.seedFinderConfig.m_layerGeometry, inputConnector);

}  // this is not Gbts config type because it is a member of the algs config,
   // which is of type Gbts cofig

// execute:
ActsExamples::ProcessCode ActsExamples::GbtsSeedingAlgorithm::execute(
    const AlgorithmContext &ctx) const {
      
      std::cout<<"Jasper: Calling the Execute member functions of class GbtsSeedingAlgorithm from Examples/Algorithms/TrakFinding/src/GbtsSeedingAlgorithm.cpp,\n This takes an input of AlgorithmContext class and returns an instance of ProcessCode"<<std::endl; 

  std::vector<Acts::Experimental::GbtsSP<SimSpacePoint>> GbtsSpacePoints =
      MakeGbtsSpacePoints(ctx, m_cfg.ActsGbtsMap);

  for (auto sp : GbtsSpacePoints) {
    const auto &links = sp.SP->sourceLinks();
    if (!links.empty()) {
      ACTS_DEBUG("Gbts space points:  Gbts_id: "
                 << sp.gbtsID << " z: " << sp.SP->z() << " r: " << sp.SP->r()
                 << " ACTS volume:  "
                 << links.front().get<IndexSourceLink>().geometryId().volume());
    }
  }
  // this is now calling on a core algorithm
  Acts::Experimental::SeedFinderGbts<SimSpacePoint> finder(
      m_cfg.seedFinderConfig, *m_gbtsGeo,
      logger().cloneWithSuffix("GbtdFinder"));

  // output of function needed for seed

  finder.loadSpacePoints(GbtsSpacePoints);

  // trigGbts file :
  Acts::Experimental::RoiDescriptor internalRoi(
      0, -4.5, 4.5, 0, -std::numbers::pi, std::numbers::pi, 0, -150., 150.);
  // ROI file:
  //  Acts::Experimental::RoiDescriptor internalRoi(0, -5, 5, 0,
  //  -std::numbers::pi, std::numbers::pi, 0, -225., 225.);

  // new version returns seeds
  SimSeedContainer seeds = finder.createSeeds(internalRoi, *m_gbtsGeo);

  m_outputSeeds(ctx, std::move(seeds));

  return ActsExamples::ProcessCode::SUCCESS;
}

std::map<std::pair<int, int>, std::pair<int, int>>
ActsExamples::GbtsSeedingAlgorithm::makeActsGbtsMap() const {

  std::cout<<"Jasper: Calling the makeActsGbtsMap member function of class GbtsSeedingAlgorithm Examples/Algorithms/TrakFinding/src/GbtsSeedingAlgorithm.cpp,\n no input and returns a map of two pairs of integers "<<std::endl; 
  std::map<std::pair<int, int>, std::pair<int, int>> ActsGbts; //define the map object
  std::ifstream data(
      m_cfg.layerMappingFile);  // read the csv file in 0 in this file refers to no Gbts ID

  std::string line; //string that will hold line of data 

  std::vector<std::vector<std::string>> parsedCsv;

  while (std::getline(data, line)) {//whilst there are still lines to read, assign the line to the string "line"

    std::stringstream lineStream(line); //reads in the string to be manipulated 
    std::string cell;
    std::vector<std::string> parsedRow; //the vector that holds all the parsed information of that line

    while (std::getline(lineStream, cell, ',')) {//while there is still info to be parsed on that line, parse each set of info seperated by a ","
      parsedRow.push_back(cell);//add to vector 
    }

    parsedCsv.push_back(parsedRow);//add vector of all the parsed info from a line to the parsed csv. indexed on line and info in line 
  }
  // file in format ACTS_vol,ACTS_lay,ACTS_mod,Gbts_id
  for (auto i : parsedCsv) {//for each vector of parsed line info

    //get the ACTS volume, layer and mod label as well as the gbts (or combined) ID and eta mod 
    int ACTS_vol = stoi(i[0]);
    int ACTS_lay = stoi(i[1]);
    int ACTS_mod = stoi(i[2]);
    int Gbts = stoi(i[5]);
    int eta_mod = stoi(i[6]);
    //create combined acts id
    int ACTS_joint = ACTS_vol * 100 + ACTS_lay;
    ActsGbts.insert({{ACTS_joint, ACTS_mod}, {Gbts, eta_mod}}); //create map of acts to gbts style combined id and mod 
  }

  return ActsGbts;
}

std::vector<Acts::Experimental::GbtsSP<ActsExamples::SimSpacePoint>>
ActsExamples::GbtsSeedingAlgorithm::MakeGbtsSpacePoints(
    const AlgorithmContext &ctx,
    std::map<std::pair<int, int>, std::pair<int, int>> map) const {

      std::cout<<"Jasper: Calling the MakeGbtsSpacePoints member function of class GbtsSeedingAlgorithm from Examples/Algorithms/TrakFinding/src/GbtsSeedingAlgorithm.cpp,\n inputs is the AlgorithmContext class and a map of pairs of ints. It returns an instance of the GbtsSP class"<<std::endl; 
  // create space point vectors
  std::vector<Acts::Experimental::GbtsSP<ActsExamples::SimSpacePoint>>
      gbtsSpacePoints;
  gbtsSpacePoints.reserve(
      m_inputSpacePoints.size());  // not sure if this is enough

  // for loop filling space
  for (const auto &isp : m_inputSpacePoints) {
    for (const auto &spacePoint : (*isp)(ctx)) {
      // Gbts space point vector
      // loop over space points, call on map
      const auto &sourceLink = spacePoint.sourceLinks();
      const auto &indexSourceLink = sourceLink.front().get<IndexSourceLink>();

      // warning if source link empty
      if (sourceLink.empty()) {
        // warning in officaial acts format
        ACTS_WARNING("warning source link vector is empty");
        continue;
      }
      int ACTS_vol_id = indexSourceLink.geometryId().volume();
      int ACTS_lay_id = indexSourceLink.geometryId().layer();
      int ACTS_mod_id = indexSourceLink.geometryId().sensitive();

      // dont want strips or HGTD
      if (ACTS_vol_id == 2 || ACTS_vol_id == 22 || ACTS_vol_id == 23 ||
          ACTS_vol_id == 24) {
        continue;
      }

      // Search for vol, lay and module=0, if doesn't exist (end) then search
      // for full thing vol*100+lay as first number in pair then 0 or mod id
      //why?
      auto ACTS_joint_id = ACTS_vol_id * 100 + ACTS_lay_id;
      auto key = std::make_pair(
          ACTS_joint_id,
          0);  // here the key needs to be pair of(vol*100+lay, 0)
      auto Find = map.find(key);

      if (Find ==
          map.end()) {  // if end then make new key of (vol*100+lay, modid)
        key = std::make_pair(ACTS_joint_id, ACTS_mod_id);  // mod ID
        Find = map.find(key);
      }

      // warning if key not in map
      if (Find == map.end()) {
        ACTS_WARNING("Key not found in Gbts map for volume id: "
                     << ACTS_vol_id << " and layer id: " << ACTS_lay_id);
        continue;
      }

      // now should be pixel with Gbts ID:
      int Gbts_id =
          Find->second
              .first;  // new map the item is a pair so want first from it

      if (Gbts_id == 0) {
        ACTS_WARNING("No assigned Gbts ID for key for volume id: "
                     << ACTS_vol_id << " and layer id: " << ACTS_lay_id);
      }

      // access IDs from map
      int eta_mod = Find->second.second;
      int combined_id = Gbts_id * 1000 + eta_mod;

      float ClusterWidth =
          0;  // false input as this is not available in examples
      // fill Gbts vector with current sapce point and ID
      gbtsSpacePoints.emplace_back(&spacePoint, Gbts_id, combined_id,
                                   ClusterWidth);  // make new GbtsSP here !
    }
  }
  ACTS_VERBOSE("Space points successfully assigned Gbts ID");

  return gbtsSpacePoints;
}

//this is what i need to compare to the athena version
std::vector<Acts::Experimental::TrigInDetSiLayer>
ActsExamples::GbtsSeedingAlgorithm::LayerNumbering() const {

  std::cout<<"Jasper: Calling the LayerNumbering member function of class GbtsSeedingAlgorithm Examples/Algorithms/TrakFinding/src/GbtsSeedingAlgorithm.cpp,\n it doesnt take inputs and returns a vector of instances of Acts::Experimental::TrigInDetSiLayer classes"<<std::endl; 
  std::vector<Acts::Experimental::TrigInDetSiLayer> input_vector; //vector will be outputted 
  
  std::vector<std::size_t> count_vector;

  m_cfg.trackingGeometry->visitSurfaces([this, &input_vector, &count_vector]( //we are visiting every surface that gives a readout in the acts version of the ITk geometry 
                                            const Acts::Surface *surface) {//for each surface (aka each actual module)
    
    Acts::GeometryIdentifier geoid = surface->geometryId();//get its acts geoemtryID's for each surface
    auto ACTS_vol_id = geoid.volume();
    auto ACTS_lay_id = geoid.layer();
    auto mod_id = geoid.sensitive();
    //actual geometry of the modules (centre and bounds)
    auto bounds_vect = surface->bounds().values();
    auto center = surface->center(Acts::GeometryContext());

    // convert bounds from local to global coordinates 
    Acts::Vector3 globalFakeMom(1, 1, 1);
    Acts::Vector2 min_bound_local =
        Acts::Vector2(bounds_vect[0], bounds_vect[1]);
    Acts::Vector2 max_bound_local =
        Acts::Vector2(bounds_vect[2], bounds_vect[3]);
    Acts::Vector3 min_bound_global = surface->localToGlobal(
        Acts::GeometryContext(), min_bound_local, globalFakeMom);
    Acts::Vector3 max_bound_global = surface->localToGlobal(
        Acts::GeometryContext(), max_bound_local, globalFakeMom);

    if (min_bound_global(0) >
        max_bound_global(0)) {  // checking that max is not accidently min
      min_bound_global.swap(max_bound_global);
    }

    float rc = 0.0;
    float minBound = 100000.0;
    float maxBound = -100000.0;
    //the point of this next part is to make the trigindetsilayer 
    // first convert to Gbts ID
    auto ACTS_joint_id = ACTS_vol_id * 100 + ACTS_lay_id;

    auto key = std::make_pair(ACTS_joint_id, 0); // 0 is because we are looking at barrel modules first which dont have a mod value

    auto Find = m_cfg.ActsGbtsMap.find(key); // find the elements that corresponds to that key

    int Gbts_id = 0;               // initialise first to avoid FLTUND later

    Gbts_id = Find->second.first;  // new map, item is pair want first

    if (Find == m_cfg.ActsGbtsMap.end()) {  // if end then make new key of (vol*100+lay, modid)

      key = std::make_pair(ACTS_joint_id, mod_id);  // mod ID means we have end cap pixels 
      Find = m_cfg.ActsGbtsMap.find(key);
      Gbts_id = Find->second.first;
    }

    short barrel_ec = 0;  // a variable that says if barrrel, 0 = barrel
    int eta_mod = Find->second.second;

    // assign barrel_ec depending on Gbts_layer ID
    if (79 < Gbts_id && Gbts_id < 85) {  // 80s, barrel
      barrel_ec = 0;
    } else if (89 < Gbts_id && Gbts_id < 99) {  // 90s positive
      barrel_ec = 2;
    } else {  // 70s negative
      barrel_ec = -2;
    }
    // if we are looking at barrel we use either r as rc and z as bounds otherise (end caps) it is otherway round 
    if (barrel_ec == 0) {
      rc = sqrt(center(0) * center(0) +
                center(1) * center(1));  // barrel center in r
      // bounds of z
      if (min_bound_global(2) < minBound) {
        minBound = min_bound_global(2);
      }
      if (max_bound_global(2) > maxBound) {
        maxBound = max_bound_global(2);
      }
    } else {
      rc = center(2);  // not barrel center in Z
      // bounds of r
      float min = sqrt(min_bound_global(0) * min_bound_global(0) +
                       min_bound_global(1) * min_bound_global(1));
      float max = sqrt(max_bound_global(0) * max_bound_global(0) +
                       max_bound_global(1) * max_bound_global(1));
      if (min < minBound) {
        minBound = min;
      }
      if (max > maxBound) {
        maxBound = max;
      }
    }

    int combined_id = Gbts_id * 1000 + eta_mod; //create combined ID
    
    //creates an iterrator pointing to the object in the vector that has the combined ID equal to that of the current module
    auto current_index = 
    find_if(input_vector.begin(), input_vector.end(), 
            [combined_id](auto n) { return n.m_subdet == combined_id; });

    //if the iterrator is not the end (therefore the Trigindetsilayer exists), we add the geometry info of the module to it 
    if (current_index != input_vector.end()) {  
      std::size_t index = std::distance(input_vector.begin(), current_index);
      input_vector[index].m_refCoord += rc;
      input_vector[index].m_minBound += minBound;
      input_vector[index].m_maxBound += maxBound;
      count_vector[index] += 1;  // increase count at the index (number of modules currently added to the trigindetsilayer)

    } else {  //if there is no trigindetsilayer that corresponds to the combined ID. make one 
      //this is where i think ill need to add code in to get the "layer number" variable implemented 
      Acts::Experimental::TrigInDetSiLayer new_Gbts_ID(combined_id, barrel_ec,
                                                       rc, minBound, maxBound);

      //as both are pushed back ino the vector, they both should have the same index                                                
      input_vector.push_back(new_Gbts_ID);
      count_vector.push_back(1);  // so the element exists and not divinding by 0
    }
    //just tells us the info for each of the acts ID'd detector modules 
    if (m_cfg.fill_module_csv) {
      std::fstream fout;
      fout.open("ACTS_modules.csv",
                std::ios::out | std::ios::app);  // add to file each time
      // print to csv for each module, no repeats so dont need to make
      // map for averaging
      fout << ACTS_vol_id << ", "                                  // vol
           << ACTS_lay_id << ", "                                  // lay
           << mod_id << ", "                                       // module
           << Gbts_id << ","                                       // Gbts id
           << eta_mod << ","                                       // eta_mod
           << center(2) << ", "                                    // z
           << sqrt(center(0) * center(0) + center(1) * center(1))  // r
           << "\n";
    }
  });

  for (std::size_t i = 0; i < input_vector.size(); i++) {// for every trigindetsilayer within the vector
    // take the cumulative rc cooridnate currently inputted and divide it by how many actual modules are in that logical layer (share the same combined ID)
    input_vector[i].m_refCoord = input_vector[i].m_refCoord / count_vector[i];
  }

  return input_vector; //return the completed vector of trigindetsilayers :)
}
