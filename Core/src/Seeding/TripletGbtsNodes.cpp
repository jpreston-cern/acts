// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/TripletGbtsNodes.hpp"

#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/SpacePointFormation/detail/StripSpacePointCalibrationImpl.hpp"

namespace Acts::Experimental {

void TripletGbtsNodes::reserve(const std::size_t nNodes,
                               const std::size_t nStrips) {
  m_spacePointIndices.reserve(nNodes);
  m_positions.reserve(nNodes);
  m_layers.reserve(nNodes);
  m_variances.reserve(nNodes);
  m_stripIndex.reserve(nNodes);
  m_rawStrips.reserve(nStrips);
  m_strips.reserve(nStrips);
}

std::uint32_t TripletGbtsNodes::insert(
    const SpacePointIndex index, const float x, const float y, const float z,
    const float r, const GbtsLayerIndex layer,
    const OuterStripSpacePointCalibrationDetails* strip, const float varianceR,
    const float varianceZ) {
  const auto node = static_cast<std::uint32_t>(m_positions.size());
  m_spacePointIndices.push_back(index);
  m_positions.push_back({x, y, z, r});
  m_layers.push_back(layer);
  m_variances.push_back({varianceR, varianceZ});
  if (strip != nullptr) {
    m_stripIndex.push_back(static_cast<std::uint32_t>(m_rawStrips.size()));
    m_rawStrips.push_back(*strip);
    // derived once here rather than once per pair in the filter
    m_strips.push_back(
        Acts::detail::deriveOuterStripSpacePointCalibrationDetails(*strip));
  } else {
    m_stripIndex.push_back(detail::kNoStrip);
  }
  return node;
}

void TripletGbtsNodes::extend(
    const SpacePointContainer& spacePoints,
    const ConstSpacePointColumnProxy<GbtsLayerIndex>& layerColumn,
    const GbtsGeometry& geometry) {
  const bool strips =
      spacePoints.hasColumns(SpacePointColumns::StripCalibrationDetails);
  const bool variances = spacePoints.hasColumns(SpacePointColumns::VarianceR |
                                                SpacePointColumns::VarianceZ);
  std::size_t nStrips = 0;
  if (strips) {
    for (const auto& sp : spacePoints) {
      nStrips += geometry.layerDescription(sp.extra(layerColumn)).technology ==
                         GbtsLayerTechnology::Strip
                     ? 1
                     : 0;
    }
  }
  reserve(size() + spacePoints.size(), m_rawStrips.size() + nStrips);
  for (const auto& sp : spacePoints) {
    const GbtsLayerIndex layer = sp.extra(layerColumn);
    const bool strip = strips && geometry.layerDescription(layer).technology ==
                                     GbtsLayerTechnology::Strip;
    insert(sp.index(), sp.x(), sp.y(), sp.z(), sp.r(), layer,
           strip ? &sp.outerStripCalibrationDetails() : nullptr,
           variances ? sp.varianceR() : 0.f, variances ? sp.varianceZ() : 0.f);
  }
}

}  // namespace Acts::Experimental
