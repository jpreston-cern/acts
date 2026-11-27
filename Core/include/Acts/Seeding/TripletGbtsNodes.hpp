// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/EventData/StripSpacePointCalibrationDetails.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace Acts::Experimental {

class GbtsGeometry;

/// The nodes of the triplet GBTS graph (`TripletGbtsGraphBuilder`): the space
/// points in the order they were inserted, node n being the n-th one, with
/// what the grid triplet search, the pixel extension and the tracking filter
/// read of them.
///
/// Unlike `GbtsNodeStorage` it keeps no eta bins, no phi order and none of the
/// per-node data of the prompt graph: the grid triplet search bins the strip
/// nodes itself, and the pixel extension the pixel nodes.
class TripletGbtsNodes final {
 public:
  /// Reserve the space for this many nodes, of which this many carry a stereo
  /// pair, so that inserting them one by one does not reallocate.
  /// @param nNodes Number of nodes to be inserted at most
  /// @param nStrips Number of them with a stereo pair at most
  void reserve(std::size_t nNodes, std::size_t nStrips);

  /// Add a node.
  /// @param index Space point index the seeds are made of
  /// @param x Node x position
  /// @param y Node y position
  /// @param z Node z position
  /// @param r Node radius, hypot(x, y)
  /// @param layer Dense GBTS layer index
  /// @param strip Stereo pair of a strip space point, null for a pixel one
  /// @param varianceR Variance of the space point in r
  /// @param varianceZ Variance of the space point in z
  /// @return The node's index
  std::uint32_t insert(SpacePointIndex index, float x, float y, float z,
                       float r, GbtsLayerIndex layer,
                       const OuterStripSpacePointCalibrationDetails* strip,
                       float varianceR, float varianceZ);

  /// Insert every space point of a container, in container order. A
  /// container carrying `SpacePointColumns::StripCalibrationDetails` has the
  /// stereo pairs of the nodes on the geometry's strip layers taken with it,
  /// and one carrying `VarianceR` and `VarianceZ` the variances (zero
  /// otherwise).
  /// @param spacePoints The space points to insert
  /// @param layerColumn Column holding the dense GBTS layer index
  /// @param geometry The GBTS geometry, which tells the strip layers
  void extend(const SpacePointContainer& spacePoints,
              const ConstSpacePointColumnProxy<GbtsLayerIndex>& layerColumn,
              const GbtsGeometry& geometry);

  /// @return Number of nodes
  std::uint32_t size() const {
    return static_cast<std::uint32_t>(m_positions.size());
  }

  /// @return The space point index of each node
  std::span<const SpacePointIndex> spacePointIndices() const {
    return m_spacePointIndices;
  }

  /// View of the node positions, layers and stereo pairs, as the tracking
  /// filter and the seeder read them. Its strip index is empty when no node
  /// carries a stereo pair.
  /// @return The view
  detail::GbtsNodeView nodeView() const {
    return detail::GbtsNodeView{
        m_positions, m_layers, m_strips,
        m_strips.empty() ? std::span<const std::uint32_t>() : m_stripIndex};
  }

  /// @return The stereo pairs as inserted, reached through `stripIndices`
  std::span<const OuterStripSpacePointCalibrationDetails> rawStrips() const {
    return m_rawStrips;
  }

  /// @return Per node the index into `rawStrips`, `detail::kNoStrip` for a
  ///         pixel node
  std::span<const std::uint32_t> stripIndices() const { return m_stripIndex; }

  /// @return The variances in r and z of each node
  std::span<const std::array<float, 2>> variances() const {
    return m_variances;
  }

 private:
  std::vector<SpacePointIndex> m_spacePointIndices;
  /// (x, y, z, r) per node
  std::vector<std::array<float, 4>> m_positions;
  std::vector<GbtsLayerIndex> m_layers;
  std::vector<std::array<float, 2>> m_variances;
  /// The stereo pairs as inserted, which the grid triplet search reads
  std::vector<OuterStripSpacePointCalibrationDetails> m_rawStrips;
  /// The same, derived once, which the tracking filter reads
  std::vector<OuterStripSpacePointCalibrationDetailsDerived> m_strips;
  /// Per node the index into both, `detail::kNoStrip` for a pixel node
  std::vector<std::uint32_t> m_stripIndex;
};

}  // namespace Acts::Experimental
