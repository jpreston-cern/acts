// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsGridTripletSearch.hpp"
#include "Acts/Seeding/TripletGbtsNodes.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"
#include "Acts/Seeding/detail/TripletGbtsGraphTypes.hpp"

#include <array>
#include <cstdint>
#include <span>
#include <tuple>
#include <vector>

namespace Acts::Experimental {

/// Builds the large radius GBTS graph from the triplets of the grid triplet
/// seeder on the strip nodes, in place of the GBTS doublet and triplet search.
///
/// A graph link is a triplet: the edge from the middle to the outer space
/// point has the edge from the inner to the middle space point as an inner
/// neighbour. Each distinct pair of space points in the triplets becomes one
/// edge. The graph can then be extended from its strip edges into the pixel
/// nodes. `GraphBasedTrackSeeder::createSeedsFromGraph` then follows the links
/// into chains, fits them with the tracking filter and removes the clones, as
/// for a graph GBTS built itself.
///
/// The nodes are a `TripletGbtsNodes`, not the prompt graph's
/// `GbtsNodeStorage`. The builder is cheap and may report per-event counts:
/// make one per event.
class TripletGbtsGraphBuilder {
 public:
  using EdgeType = detail::TripletGbtsEdge;

  struct Config {
    /// Most inner neighbours an edge links to, the best triplets first; 0 for
    /// no limit
    std::uint32_t maxInnerNeighbours = detail::kTripletGbtsMaxEdgeNeighbours;

    /// Extend the graph into the pixel nodes (those without a stereo pair):
    /// from each strip edge that starts a triplet, inner pixel nodes whose
    /// triplet with the edge matches the edge's own triplets, and from each
    /// such pixel to strip edge, inner pixel nodes again. The triplets then
    /// need only hold strip space points.
    bool pixelExtension = false;
    /// Largest radius of the inner strip node of an edge extended into pixels
    float pixelMaxStripR = 600.f;
    /// Azimuthal window around the outer node of the new edge (rad)
    float pixelPhiWindow = 0.06f;
    /// Smallest radial distance of a pixel node to the node it links to (mm)
    float pixelMinDeltaR = 10.f;
    /// Curvature match to the edge's triplets: absolute (1/mm) plus relative
    /// to the larger curvature
    float pixelMaxDeltaCurvature = 2e-4f;
    float pixelMaxRelCurvature = 0.2f;
    /// dz/dr match of the new edge to the edge: absolute plus relative to it
    float pixelMaxDeltaTau = 0.05f;
    float pixelMaxRelDeltaTau = 0.05f;
    /// Most pixel nodes linked to per edge, the best curvature match first
    std::uint32_t pixelMaxPerEdge = 5;
    /// Edges with |eta| (from their own dz/dr) below this use the barrel
    /// matching parameters below instead; 0 for one set everywhere
    float pixelBarrelMaxAbsEta = 0.f;
    float pixelPhiWindowBarrel = 0.06f;
    float pixelMaxDeltaCurvatureBarrel = 2e-4f;
    float pixelMaxRelCurvatureBarrel = 0.2f;
    float pixelMaxDeltaTauBarrel = 0.05f;
    float pixelMaxRelDeltaTauBarrel = 0.05f;
    std::uint32_t pixelMaxPerEdgeBarrel = 5;
    /// Extend an edge only into the pixel layers that the geometry's layer
    /// connection table connects the layer of its inner node to
    bool pixelLayersFromTable = false;
    /// Keep only the grid triplets whose two node pairs are both connections
    /// of the geometry's layer connection table, so that the table decides
    /// which strip layers are linked, as it does for the pixel extension
    bool stripLayersFromTable = false;
    /// Let the grid search combine only the layers the geometry's layer
    /// connection table connects: each grid bin is split by layer, and a
    /// middle layer is searched with only its connected bottom and top layers
    /// in the grid's neighbour bins. The triplets `stripLayersFromTable`
    /// would drop are then never made, so the filter can keep others instead.
    bool stripSearchFromTable = false;
    /// Let the layer connection table decide which strip edges are extended
    /// into pixels: an edge is extended when the table connects the layer of
    /// its inner node to a pixel layer, in place of `pixelMaxStripR`. The pixel
    /// layers come from the table as with `pixelLayersFromTable`;
    /// `pixelMinDeltaR` still applies, as the slope and curvature of a pair of
    /// nodes close in radius are dominated by their measurement errors.
    bool pixelLinksFromTable = false;
  };

  /// Counts of the last graph built
  struct Stats {
    std::uint32_t triplets = 0;
    std::uint32_t edges = 0;
    std::uint32_t links = 0;
    std::uint32_t linksOverCap = 0;
    std::uint32_t edgesLeftOut = 0;
    std::uint32_t pixelStripLinks = 0;
    std::uint32_t pixelPixelLinks = 0;
    /// Grid triplets dropped by `stripLayersFromTable`
    std::uint32_t tripletsOutsideTable = 0;
    /// Grid triplets dropped by the search's forward quality cut
    std::uint32_t tripletsBelowQuality = 0;
  };

  /// @param config The configuration
  /// @param search The triplet search, which must outlive the builder
  /// @param geometry The GBTS geometry, whose layer connection table the
  ///        table options read; must outlive the builder
  /// @param offset Subtracted from the node x and y positions for the search,
  ///        e.g. the beam spot when the nodes were not corrected for it
  /// @param stats Filled with the counts of the graph built, if not null
  TripletGbtsGraphBuilder(const Config& config,
                          const GbtsGridTripletSearch& search,
                          const GbtsGeometry& geometry,
                          std::array<float, 2> offset, Stats* stats = nullptr)
      : m_cfg(config),
        m_search(&search),
        m_geometry(&geometry),
        m_offset(offset),
        m_stats(stats) {}

  /// Search the triplets among the strip nodes and build the graph from them
  /// @param nodes The nodes; the strip ones enter the search in node order
  /// @param bFieldInZ The magnetic field along z, for the triplet search
  /// @return The graph, ordered outside in
  detail::GbtsGraph<EdgeType> buildTheGraph(const TripletGbtsNodes& nodes,
                                            float bFieldInZ) const;

 private:
  template <typename edge_of_t>
  void extendIntoPixels(
      const detail::GbtsNodeView& nodeView, std::uint32_t nNodes,
      std::vector<EdgeType>& edges, edge_of_t& edgeOf,
      std::vector<std::tuple<std::uint32_t, float, std::uint32_t>>& rawLinks,
      std::span<const std::uint8_t> allowedPairs, std::uint32_t nLayers,
      Stats& stats) const;

  Config m_cfg;
  const GbtsGridTripletSearch* m_search = nullptr;
  const GbtsGeometry* m_geometry = nullptr;
  std::array<float, 2> m_offset{0.f, 0.f};
  Stats* m_stats = nullptr;
};

}  // namespace Acts::Experimental
