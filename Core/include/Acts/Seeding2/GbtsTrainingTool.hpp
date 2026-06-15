// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// includes needed
#pragma once
#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace Acts::Experimental {

struct TrackCoordinates {

  float r{};
  float z{};
};

struct LayerDescription {
  LayerDescription(float minR_, float maxR_, float minZ_, float maxZ_,
                   std::int32_t gbtsId_);

  // r values
  float minR{};
  float maxR{};

  // z values
  float minZ{};
  float maxZ{};

  std::int32_t gbtsId{};
};

using LayerIdPair = std::pair<std::int32_t, std::int32_t>;

struct LayerIdPairHash {
  std::size_t operator()(const LayerIdPair& pair) const noexcept {
    const auto h1 = std::hash<std::int32_t>{}(pair.first);
    const auto h2 = std::hash<std::int32_t>{}(pair.second);

    return h1 ^ (h2 << 1);
  }
};

class GbtsTrainingTool {
 public:
  explicit GbtsTrainingTool(std::string& geometryInformation);

  void addTrack(const std::vector<TrackCoordinates>& Track);

  void createConnectionTable(const std::string& outputFileLocations,
                             const double probThreshold) const;
  
 std::optional<std::int32_t> findGbtsIdByCoord(float r, float z) const;

 private:
  
  std::vector<LayerDescription> m_detectorGeometry{};

  std::unordered_map<LayerIdPair, std::uint32_t, LayerIdPairHash> m_layerPairs{};

  std::uint32_t m_totalTracks = 0;
};

}  // namespace Acts::Experimental
