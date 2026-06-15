// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/GbtsTrainingTool.hpp"

#include <fstream>
#include <iostream>
#include <string>

namespace Acts::Experimental {

LayerDescription::LayerDescription(float minR_, float maxR_, float minZ_,
                                   float maxZ_, std::int32_t gbtsId_)
    : minR(minR_), maxR(maxR_), minZ(minZ_), maxZ(maxZ_), gbtsId(gbtsId_) {}

GbtsTrainingTool::GbtsTrainingTool(std::string& geometryInformation) {

  std::ifstream inStream(geometryInformation.c_str());
  // define how many lines there are for reserving
  std::uint32_t lines{};
  std::string line{};
  while (std::getline(inStream, line)) {
    lines++;
  }
  inStream.clear();
  inStream.seekg(0);

  //reserves
  m_detectorGeometry.reserve(lines);
  m_layerPairs.reserve(lines * (lines - 1));
  
  // create geometry objects
  float minR{};
  float maxR{};

  float minZ{};
  float maxZ{};

  std::int32_t gbtsId{};

  for (std::uint32_t l = 0; l < lines; l++) {
    inStream >> minR >> maxR >> minZ >> maxZ >> gbtsId;

    m_detectorGeometry.emplace_back(minR, maxR, minZ, maxZ, gbtsId);
  }

  // create map linked pairs of GBTS ids with number of transitions between
  // layers
  for (std::uint32_t i = 0; i < m_detectorGeometry.size(); i++) {
    for (std::uint32_t j = 0; j < m_detectorGeometry.size(); j++) {

      if(i == j){
        continue;
      }

      LayerIdPair pair;
      pair.first = m_detectorGeometry[i].gbtsId;
      pair.second = m_detectorGeometry[j].gbtsId;

      // key = GBTS ids of pair, value = number of transitions
      m_layerPairs.emplace(pair, 0);
    }
  }
}

void GbtsTrainingTool::addTrack(const std::vector<TrackCoordinates>& track) {
  if (track.size() < 2) {
    std::cerr << "Warning: Track only has one measurement, skipping" << std::endl;
    return;
  }
  // container for gbts IDs of the track
  std::vector<std::int32_t> layerGbtsIds{};
  layerGbtsIds.reserve(track.size());

  // find GBTS ids for all measurements in a track
  for (const auto& measurement : track) {
    const float r = measurement.r;
    const float z = measurement.z;
    
    const auto gbtsId = findGbtsIdByCoord(r, z);
    if(!gbtsId){
      std::cerr << "Warning: no Gbts Layer for coordinates with r: " << r << " and z: " << z << std::endl;
      continue;
    }
    layerGbtsIds.emplace_back(gbtsId.value());

  }

  // update map with track layer transitions
  for (std::uint32_t id = 0; id + 1 < layerGbtsIds.size(); id++) {
    const std::int32_t index1 = layerGbtsIds[id];
    const std::int32_t index2 = layerGbtsIds[id + 1];

    if (index1 == index2) {
      std::cout << "Warning: track transitions between same layer, skipping" << std::endl;

      continue;
    }

    m_layerPairs[{index1, index2}] += 1;

  }

  m_totalTracks++;
  // std::ofstream hitOutput;
  // std::ofstream mapOutput;
  // if (m_totalTracks == 1){

  //   hitOutput.open("hit_output.txt");
  //   mapOutput.open("map_output.txt");
  //   for(const auto& hit : track){
  //     hitOutput << "r: " << hit.r << " z: " << hit.z << "\n";
  //   }
  //   for(const auto& [pair, nTransitions] : m_layerPairs){
  //     mapOutput << "Layer Pair: " << pair.first << ", " << pair.second << " Transitions: " << nTransitions << "\n";
  //   }
  // }
}

void GbtsTrainingTool::createConnectionTable(
    const std::string& outputFileLocation,
    const double probThreshold) const {
  if (m_totalTracks == 0) {
    std::runtime_error(
        "Warning: no tracks were added when creating connection table");
    return;
  }
  // define output text file
  std::ofstream outputFile(outputFileLocation);

  // std::cout<<"Jasper: the number of  of 80000 -> 82000 transitions is: " << m_layerPairs.at({80000,82000}) << std::endl;
  // std::cout<<"Jasper the number of  of 80000 -> 82000 transitions is: " << m_layerPairs.at({80000,82000})/static_cast<float>(m_totalTracks) << std::endl;
  // std::cout<<"Jasper: total tracks are: " << m_totalTracks << std::endl; 
  // for every layer transition, only output those that pass probability cut
  for (const auto& [indexes, nTransitions] : m_layerPairs) {
    const float probability =
        static_cast<float>(nTransitions) / static_cast<float>(m_totalTracks);

    if (probThreshold < probability) {
      outputFile << indexes.first << " " << indexes.second << "\n";
    }
  }
}

std::optional<std::int32_t> GbtsTrainingTool::findGbtsIdByCoord(const float r,
                                                 const float z) const {

  for (const auto& layer : m_detectorGeometry) {

    const float zTolMin = layer.minZ - 2.0f;
    const float zTolMax = layer.maxZ + 2.0f;
    const float rTolMin = layer.minR - 2.0f;
    const float rTolMax = layer.maxR + 2.0f;

    if (zTolMin <= z && z <= zTolMax) {
      if (rTolMin <= r && r <= rTolMax) {
        return layer.gbtsId;
      }
    }
  }

  return std::nullopt;
}

}  // namespace Acts::Experimental
