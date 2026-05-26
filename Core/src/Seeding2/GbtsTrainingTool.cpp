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

GbtsTrainingTool::GbtsTrainingTool(std::istream& inStream) {
  // define how many lines there are for reserving
  std::uint32_t lines{};
  std::string line{};
  while (std::getline(inStream, line)) {
    lines++;
  }
  inStream.clear();
  inStream.seekg(0);

  m_detectorGeometry.reserve(lines);

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
    for (std::uint32_t j = i + 1; j < m_detectorGeometry.size(); j++) {
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
    std::cout << "Warning: Track only has one measurement" << std::endl;
    return;
  }
  // container for gbts IDs of the track
  std::vector<std::int32_t> layerGbtsIds{};
  layerGbtsIds.reserve(track.size());

  // find GBTS ids for all measurements in a track
  for (const auto& measurement : track) {
    const float r = measurement.r;
    const float z = measurement.z;

    const std::int32_t gbtsId = findGbtsIdByCoord(r, z);

    layerGbtsIds.emplace_back(gbtsId);
  }

  // update map with track layer transitions
  for (std::uint32_t id = 0; id + 1 < layerGbtsIds.size(); id++) {
    const std::int32_t index1 = layerGbtsIds[id];
    const std::int32_t index2 = layerGbtsIds[id + 1];

    if (index1 == index2) {
      std::cout << "Warning: track transitions between same layer" << std::endl;

      continue;
    }

    m_layerPairs[{index1, index2}] += 1;
  }

  m_totalTracks += 1;
}

void GbtsTrainingTool::createConnectionTable(
    const std::filesystem::path& outputFileLocation,
    const double probThreshold) const {
  if (m_totalTracks == 0) {
    std::runtime_error(
        "Warning: no tracks were added when creating connection table");
    return;
  }
  // define output text file
  std::ofstream outputFile(outputFileLocation);

  // for every layer transition, only output those that pass probability cut
  for (const auto& [indexes, nTransitions] : m_layerPairs) {
    const double probability =
        static_cast<double>(nTransitions) / static_cast<double>(m_totalTracks);

    if (probThreshold < probability) {
      outputFile << indexes.first << " " << indexes.second << "\n";
    }
  }
}

std::int32_t GbtsTrainingTool::findGbtsIdByCoord(const float r,
                                                 const float z) const {
  for (const auto& layer : m_detectorGeometry) {
    if (layer.minZ < z && z < layer.maxZ) {
      if (layer.minR < r && r < layer.maxR) {
        return layer.gbtsId;
      }
    }
  }

  throw std::runtime_error("No matching GBTS ID found for coordinates");
}

}  // namespace Acts::Experimental
