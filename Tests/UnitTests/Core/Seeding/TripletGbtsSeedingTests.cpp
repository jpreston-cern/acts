// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding/DoubletSeedFinder.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsGridTripletSearch.hpp"
#include "Acts/Seeding/GbtsLayerConnection.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"
#include "Acts/Seeding/TripletGbtsGraphBuilder.hpp"
#include "Acts/Seeding/TripletGbtsNodes.hpp"
#include "Acts/Seeding/TripletGbtsSeedSelector.hpp"
#include "Acts/Seeding/TripletSeedFinder.hpp"
#include "Acts/Seeding/detail/TripletGbtsGraphTypes.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Harness for the triplet GBTS graph, the displaced (large radius) graph built
// from grid triplets: displaced helices on a toy strip barrel, a grid triplet
// search on the graph nodes, and checks, from the geometry alone, that the
// graph and the seeds extracted from it hold every track.

namespace Acts::Test {

namespace {

using namespace Acts::UnitLiterals;
using Experimental::GbtsLayerTechnology;
using Experimental::GbtsLayerType;

/// The region the triplet graph is asked to cover.
constexpr float kD0Max = 300.f;
constexpr float kZ0Max = 500.f;
constexpr float kMaxOuterRadius = 1050.f;

/// Half-length in z of every strip barrel layer. Wide enough that a track at
/// the z0 limit still crosses all four.
constexpr float kBarrelHalfZ = 1500.f;

/// Eta bin width the layers are split into, as in ATLAS'
/// createLinkingScheme.py.
constexpr float kEtaBinWidth = 0.2f;

/// One strip barrel layer. The ids follow ATLAS, which the seeder still keys
/// on: barrel ids 1000 apart are adjacent.
struct LayerSpec {
  Experimental::GbtsExperimentLayerId id{};
  float radius{};
};

/// Layers plus the connector links between them, outer (src) to inner (dst).
struct ToyDetector {
  std::vector<LayerSpec> layers;
  std::vector<std::pair<Experimental::GbtsExperimentLayerId,
                        Experimental::GbtsExperimentLayerId>>
      links;
  /// Give the layers a depth, innermost first, for the innermost-layer cuts.
  /// Off, as the experiment leaves it for its displaced seeding.
  bool withDepth = false;
};

/// Four strip barrel layers, spaced roughly like the ITk strip barrel, all
/// inside `kMaxOuterRadius`.
///
/// Every layer sits outside `kD0Max`, which the displaced phi window needs:
/// it takes acos(d0Max / r) of a bin's radius, and a layer inside the d0
/// limit would hand that a value above one.
ToyDetector stripBarrelDetector() {
  return {{{80000, 400.f}, {81000, 560.f}, {82000, 760.f}, {83000, 1000.f}},
          {{81000, 80000}, {82000, 81000}, {83000, 82000}}};
}

/// Changes the layer pairs after they are filled.
using ConnectionTweak =
    std::function<void(std::vector<Experimental::GbtsLayerConnection>&)>;

std::shared_ptr<Experimental::GbtsGeometry> makeGeometry(
    const ToyDetector& detector, const ConnectionTweak& connectionTweak = {}) {
  std::vector<Experimental::GbtsLayerDescription> layers;
  layers.reserve(detector.layers.size());
  for (const LayerSpec& spec : detector.layers) {
    Experimental::GbtsLayerDescription layer;
    layer.id = spec.id;
    layer.type = GbtsLayerType::Barrel;
    layer.technology = GbtsLayerTechnology::Strip;
    layer.refCoord = spec.radius;
    layer.minBound = -kBarrelHalfZ;
    layer.maxBound = kBarrelHalfZ;
    if (detector.withDepth) {
      layer.depth = static_cast<std::int32_t>(layers.size());
    }
    layers.push_back(layer);
  }

  std::vector<Experimental::GbtsLayerConnection> connections;
  connections.reserve(detector.links.size());
  for (const auto& [src, dst] : detector.links) {
    connections.push_back({src, dst});
  }
  if (connectionTweak) {
    connectionTweak(connections);
  }

  // The bin table is what decides which eta bins may be connected at all, and
  // it is built against a z0 range. Left at its default of +-168mm it would
  // only link bins a track from the beamspot could cross, and a displaced
  // track starting half a metre down the beamline would have no link to
  // travel along however well it passed the cuts.
  const Experimental::GbtsZ0Range z0Range{-kZ0Max, kZ0Max};

  return std::make_shared<Experimental::GbtsGeometry>(layers, connections,
                                                      kEtaBinWidth, z0Range);
}

/// A displaced helix, given at its perigee.
///
/// The prompt suite's tracks are straight lines through the origin, which is
/// the one thing the triplet graph may not assume, so these carry a real
/// impact parameter and a real curvature instead.
struct Track {
  /// Signed transverse impact parameter.
  float d0{};
  /// Longitudinal impact parameter.
  float z0{};
  /// Azimuth of the track direction at the perigee.
  float phi0{};
  /// cot(theta), constant along a helix.
  float tau{};
  /// Transverse radius of curvature.
  float radius{};
  /// +1 for a track turning anticlockwise.
  int sense{};

  /// Centre of the transverse circle.
  std::array<float, 2> centre() const {
    // the left normal of the direction at the perigee, which both the perigee
    // and the centre lie along
    const float nx = -std::sin(phi0);
    const float ny = std::cos(phi0);
    const float d = d0 + static_cast<float>(sense) * radius;
    return {d * nx, d * ny};
  }

  /// Transverse position at turn angle @p t from the perigee.
  std::array<float, 2> position(const float t) const {
    const std::array<float, 2> c = centre();
    const float vx = d0 * -std::sin(phi0) - c[0];
    const float vy = d0 * std::cos(phi0) - c[1];
    const float a = static_cast<float>(sense) * t;
    return {c[0] + vx * std::cos(a) - vy * std::sin(a),
            c[1] + vx * std::sin(a) + vy * std::cos(a)};
  }

  /// Unit direction at turn angle @p t, transverse part then tau.
  std::array<float, 3> direction(const float t) const {
    const std::array<float, 2> c = centre();
    const float vx = d0 * -std::sin(phi0) - c[0];
    const float vy = d0 * std::cos(phi0) - c[1];
    const float a =
        static_cast<float>(sense) * t + 0.5f * std::numbers::pi_v<float>;
    const float dx =
        static_cast<float>(sense) * (vx * std::cos(a) - vy * std::sin(a));
    const float dy =
        static_cast<float>(sense) * (vx * std::sin(a) + vy * std::cos(a));
    return {dx / radius, dy / radius, tau};
  }

  /// z at turn angle @p t. The transverse path from the perigee is radius * t.
  float z(const float t) const { return z0 + tau * radius * t; }
};

/// Where a track crosses a layer, or nothing if it never reaches it or leaves
/// the barrel first. The transverse radius grows monotonically with the turn
/// angle over half a turn, so a bisection settles it.
std::optional<float> crossingAngle(const Track& track, const float radius) {
  const auto radiusAt = [&](const float t) {
    const std::array<float, 2> p = track.position(t);
    return std::hypot(p[0], p[1]);
  };

  constexpr float pi = std::numbers::pi_v<float>;
  if (radiusAt(pi) < radius) {
    // the track turns back before reaching the layer
    return std::nullopt;
  }

  float lo = 0.f;
  float hi = pi;
  for (int i = 0; i < 80; ++i) {
    const float mid = 0.5f * (lo + hi);
    (radiusAt(mid) < radius ? lo : hi) = mid;
  }
  return 0.5f * (lo + hi);
}

/// Tracks spread over the full region the graph is configured for: the whole
/// d0 and z0 range, both charges, both signs of eta and all of phi.
std::vector<Track> makeDisplacedTracks() {
  constexpr float pi = std::numbers::pi_v<float>;
  constexpr std::size_t nTracks = 16;

  std::vector<Track> tracks;
  tracks.reserve(nTracks);

  for (std::size_t i = 0; i < nTracks; ++i) {
    const float frac = static_cast<float>(i) / nTracks;
    Track track;
    // half a step offset so no track sits on the phi wrap-around
    track.phi0 = -pi + 2 * pi * (frac + 0.5f / nTracks);
    // the full d0 range, alternating sign, short of the configured limit
    track.d0 = (i % 2 == 0 ? 1.f : -1.f) * 0.92f * kD0Max * frac;
    // the full z0 range, walked the other way so it does not track d0
    track.z0 = (i % 4 < 2 ? 1.f : -1.f) * 0.95f * kZ0Max * (1.f - frac);
    track.tau = -0.8f + 1.6f * frac;
    // 2 to 10 GeV at 2T, comfortably inside the curvature cut
    track.radius = 3333.f + 13333.f * frac;
    track.sense = (i % 2 == 0) ? 1 : -1;
    tracks.push_back(track);
  }
  return tracks;
}

/// A tight bundle of displaced tracks, close enough in every coordinate that
/// the graph forms edges across them.
///
/// The prompt suite's dense case packs tracks into a narrow phi and tau
/// window. Displaced, phi is not enough on its own: d0 moves a track's azimuth
/// at a given radius by acos(d0 / r), so a wide d0 range fans the bundle back
/// out however close its directions are. These share a d0 to within a couple
/// of centimetres, which is what puts their hits on top of one another, and it
/// is the case the prompt graph never has to face.
///
/// The tau range is the one thing kept open. Two tracks closer in tau than
/// `tauRatioCut` are not two tracks as far as the graph is concerned, so
/// packing them tighter than that measures nothing but the fixture: adjacent
/// tracks here differ by about 1.5%, against a cut of 0.7%.
std::vector<Track> makeDenseDisplacedTracks() {
  constexpr std::size_t nTracks = 40;

  std::vector<Track> tracks;
  tracks.reserve(nTracks);

  for (std::size_t i = 0; i < nTracks; ++i) {
    const float frac = static_cast<float>(i) / nTracks;
    Track track;
    track.phi0 = -0.03f + 0.06f * frac;
    track.d0 = 240.f + 24.f * frac;
    track.z0 = 240.f + 40.f * frac;
    track.tau = 0.20f + 0.30f * frac;
    track.radius = 3333.f + 13333.f * frac;
    track.sense = (i % 2 == 0) ? 1 : -1;
    tracks.push_back(track);
  }
  return tracks;
}

/// A strip module of the toy barrel, roughly the ITk strip barrel.
constexpr float kStereoAngle = 26e-3f;
constexpr float kModuleGap = 2.f;
constexpr float kStripHalfLength = 30.f;
/// Position variances of the toy space points: about the z resolution of a
/// stereo strip pair, and a small radial one.
constexpr float kVarianceZ = 0.1f;
constexpr float kVarianceR = 0.01f;

/// One hit per crossed layer per track.
///
/// Every hit carries the stereo pair it was formed from, built so that the
/// calibration run against the true direction returns the true crossing: the
/// outer sensor is centred on it and the inner one sits a module gap back
/// along the track.
///
/// @param detector the toy detector
/// @param tracks the tracks crossing it
/// @param walk how far along the strip the uncalibrated point is put
/// @return the space points
SpacePointContainer makeSpacePoints(const ToyDetector& detector,
                                    const std::vector<Track>& tracks,
                                    const float walk = 0.f) {
  SpacePointContainer container(
      SpacePointColumns::CopiedFromIndex | SpacePointColumns::X |
      SpacePointColumns::Y | SpacePointColumns::Z | SpacePointColumns::R |
      SpacePointColumns::Phi | SpacePointColumns::VarianceZ |
      SpacePointColumns::VarianceR | SpacePointColumns::StripCalibrationDetails);

  auto layerColumn =
      container.createColumn<Experimental::GbtsLayerIndex>("gbtsLayerIndex");
  auto clusterWidthColumn = container.createColumn<float>("clusterWidth");
  auto localPositionColumn = container.createColumn<float>("localPositionY");
  auto trackColumn = container.createColumn<std::uint32_t>("trackId");

  container.reserve(
      static_cast<std::uint32_t>(tracks.size() * detector.layers.size()));

  for (std::size_t track = 0; track < tracks.size(); ++track) {
    for (std::size_t layer = 0; layer < detector.layers.size(); ++layer) {
      const std::optional<float> turn =
          crossingAngle(tracks[track], detector.layers[layer].radius);
      if (!turn.has_value()) {
        continue;
      }

      const std::array<float, 2> xy = tracks[track].position(*turn);
      const float z = tracks[track].z(*turn);
      if (std::abs(z) > kBarrelHalfZ) {
        continue;
      }
      const std::array<float, 3> point{xy[0], xy[1], z};
      const std::array<float, 3> direction = tracks[track].direction(*turn);

      // strips along z, the two sensors each rotated half the stereo angle
      // about the module normal
      const float r = std::hypot(xy[0], xy[1]);
      const std::array<float, 3> across{xy[1] / r, -xy[0] / r, 0.f};
      const float half = 0.5f * kStereoAngle;

      auto sp = container.createSpacePoint();
      // the uncalibrated point has walked along its own strip
      sp.x() = point[0] + walk * std::sin(half) * across[0];
      sp.y() = point[1] + walk * std::sin(half) * across[1];
      sp.z() = point[2] + walk * std::cos(half);
      sp.r() = std::hypot(sp.x(), sp.y());
      sp.phi() = std::atan2(sp.y(), sp.x());
      sp.copiedFromIndex() = sp.index();
      sp.varianceZ() = kVarianceZ;
      sp.varianceR() = kVarianceR;
      // the dense layer index, not the GBTS layer id
      sp.extra(layerColumn) = static_cast<Experimental::GbtsLayerIndex>(layer);
      sp.extra(clusterWidthColumn) = 0.f;
      sp.extra(localPositionColumn) = 0.f;
      sp.extra(trackColumn) = static_cast<std::uint32_t>(track);

      OuterStripSpacePointCalibrationDetails details{};
      for (std::size_t side = 0; side < 2; ++side) {
        const float sign = side == 0 ? -1.f : 1.f;
        const std::array<float, 3> axis{sign * std::sin(half) * across[0],
                                        sign * std::sin(half) * across[1],
                                        std::cos(half)};
        std::array<float, 3>& halfVector =
            side == 0 ? details.innerHalfVector : details.outerHalfVector;
        for (std::size_t i = 0; i < 3; ++i) {
          halfVector[i] = axis[i] * kStripHalfLength;
        }
      }
      for (std::size_t i = 0; i < 3; ++i) {
        // the true crossing sits at the centre of the outer strip, and the
        // inner sensor a module gap back along the track
        details.outerCenter[i] = point[i];
        details.innerToOuterSeparation[i] = kModuleGap * direction[i];
      }
      sp.outerStripCalibrationDetails() = details;
    }
  }

  return container;
}

/// The grid triplet search configured for the toy barrel and the displaced
/// region: doublets up to the outermost layer, the strip calibration on, and a
/// loose seed filter, so that the triplets sharing a doublet survive to be
/// linked.
Experimental::GbtsGridTripletSearch::Config tripletSearchConfig() {
  Experimental::GbtsGridTripletSearch::Config config;
  config.grid.minPt = 1_GeV;
  config.grid.rMin = 0.f;
  config.grid.rMax = kMaxOuterRadius + 50.f;
  config.grid.zMin = -kBarrelHalfZ - 100.f;
  config.grid.zMax = kBarrelHalfZ + 100.f;
  config.grid.deltaRMax = 700.f;
  config.grid.cotThetaMax = 5.f;
  config.grid.impactMax = kD0Max;
  config.grid.zBinEdges = {config.grid.zMin, config.grid.zMax};

  for (DoubletSeedFinder::Config* doublets :
       {&config.bottomDoublets, &config.topDoublets}) {
    doublets->spacePointsSortedByRadius = true;
    doublets->deltaRMin = 20.f;
    doublets->deltaRMax = 700.f;
    doublets->impactMax = kD0Max;
    doublets->interactionPointCut = true;
    doublets->collisionRegionMin = -kZ0Max;
    doublets->collisionRegionMax = kZ0Max;
    doublets->cotThetaMax = 5.f;
    doublets->minPt = 1_GeV;
    doublets->stripCalibrationTolerance = 1.75f;
  }
  config.bottomDoublets.candidateDirection = Direction::Backward();
  config.topDoublets.candidateDirection = Direction::Forward();

  config.triplets.useStripInfo = true;
  config.triplets.minPt = 1_GeV;
  config.triplets.impactMax = kD0Max;

  config.filter.deltaRMin = 20.f;
  // the dense bundle puts all of its tracks through each middle space point;
  // the filter's cap per middle space point is not what is tested here
  config.filter.maxSeedsPerSpM = 100;

  config.useStripDetails = true;
  return config;
}

/// The seeder and everything the triplet graph needs. The geometry, the search
/// and the statistics are held by pointer, as the graph builder keeps their
/// address.
struct SeederSetup {
  std::shared_ptr<Experimental::GbtsGeometry> geometry;
  Experimental::GraphBasedTrackSeeder seeder;
  Experimental::TripletGbtsSeedSelector selector;
  std::unique_ptr<Experimental::GbtsGridTripletSearch> search;
  std::unique_ptr<Experimental::TripletGbtsGraphBuilder::Stats> stats;
  std::unique_ptr<Experimental::TripletGbtsGraphBuilder> graph;
  Experimental::GbtsTrackingFilter filter;
  float bFieldInZ{};
};

SeederSetup makeSeeder(
    const ToyDetector& detector,
    const std::function<void(Experimental::TripletGbtsGraphBuilder::Config&)>&
        graphTweak = {},
    const std::function<void(Experimental::GraphBasedTrackSeeder::Config&,
                             Experimental::TripletGbtsSeedSelector::Config&)>&
        seederTweak = {},
    const ConnectionTweak& connectionTweak = {}) {
  auto geometry = makeGeometry(detector, connectionTweak);

  const auto makeLogger = []() -> std::unique_ptr<const Logger> {
    return getDefaultLogger("TripletGbtsTest", Logging::Level::WARNING);
  };

  Experimental::GraphBasedTrackSeeder::Config config;
  // the toy setup has no tau lookup table and no cluster widths
  config.useClusterWidthCuts = false;
  // a link here is a grid triplet, so three nodes are a seed
  config.minSeedLevel = 2;
  Experimental::TripletGbtsSeedSelector::Config selectorConfig;
  if (seederTweak) {
    seederTweak(config, selectorConfig);
  }

  Experimental::TripletGbtsGraphBuilder::Config graphConfig;
  if (graphTweak) {
    graphTweak(graphConfig);
  }

  auto search = std::make_unique<Experimental::GbtsGridTripletSearch>(
      tripletSearchConfig(), makeLogger());
  auto stats = std::make_unique<Experimental::TripletGbtsGraphBuilder::Stats>();
  auto graph = std::make_unique<Experimental::TripletGbtsGraphBuilder>(
      graphConfig, *search, *geometry, std::array<float, 2>{0.f, 0.f},
      stats.get());

  // the filter takes its z0 window from the geometry
  Experimental::GbtsTrackingFilter::Config filterConfig;

  return SeederSetup{
      .geometry = geometry,
      .seeder = Experimental::GraphBasedTrackSeeder(
          Experimental::GraphBasedTrackSeeder::DerivedConfig(config), geometry,
          makeLogger()),
      .selector =
          Experimental::TripletGbtsSeedSelector(selectorConfig, makeLogger()),
      .search = std::move(search),
      .stats = std::move(stats),
      .graph = std::move(graph),
      .filter = Experimental::GbtsTrackingFilter(filterConfig, geometry,
                                                 makeLogger()),
      .bFieldInZ = 2_T,
  };
}

/// The graph nodes of the space points.
Experimental::TripletGbtsNodes makeNodes(
    const SeederSetup& setup, const SpacePointContainer& spacePoints) {
  Experimental::TripletGbtsNodes nodes;
  nodes.extend(
      spacePoints,
      spacePoints.column<Experimental::GbtsLayerIndex>("gbtsLayerIndex"),
      *setup.geometry);
  return nodes;
}

/// The graph, plus what is needed to read it back against the tracks.
struct GraphRun {
  std::vector<Experimental::TripletGbtsGraphBuilder::EdgeType> edges;
  std::uint32_t nEdges{};
  std::uint32_t nConnections{};
  /// The inner neighbours of edge `e` are `neighbours[neighbourOffsets[e]]` up
  /// to `neighbours[neighbourOffsets[e + 1]]`, as the graph holds them.
  std::vector<std::uint32_t> neighbourOffsets;
  std::vector<std::uint32_t> neighbours;
  /// Node index of each space point, by space point index.
  std::vector<std::uint32_t> nodeOfSpacePoint;
  /// Space point index of each node.
  std::vector<std::uint32_t> spacePointOfNode;
  /// Edge index by (inner node, outer node).
  std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> edgeIndex;
  /// The builder's counts.
  Experimental::TripletGbtsGraphBuilder::Stats stats;
};

GraphRun buildGraph(const SeederSetup& setup,
                    const SpacePointContainer& spacePoints) {
  const Experimental::TripletGbtsNodes nodes = makeNodes(setup, spacePoints);

  GraphRun run;
  Experimental::detail::GbtsGraph<Experimental::TripletGbtsGraphBuilder::EdgeType>
      graph = setup.graph->buildTheGraph(nodes, setup.bFieldInZ);
  run.stats = *setup.stats;
  run.edges = std::move(graph.edgeStorage);
  run.nEdges = static_cast<std::uint32_t>(run.edges.size());
  run.nConnections = static_cast<std::uint32_t>(graph.neighbours.size());
  run.neighbourOffsets = std::move(graph.neighbourOffsets);
  run.neighbours = std::move(graph.neighbours);

  run.nodeOfSpacePoint.assign(spacePoints.size(), 0);
  run.spacePointOfNode.assign(nodes.size(), 0);
  for (std::uint32_t node = 0; node < nodes.size(); ++node) {
    run.nodeOfSpacePoint.at(nodes.spacePointIndices()[node]) = node;
    run.spacePointOfNode.at(node) = nodes.spacePointIndices()[node];
  }

  for (std::uint32_t e = 0; e < run.nEdges; ++e) {
    run.edgeIndex.emplace(std::pair{run.edges[e].n1, run.edges[e].n2}, e);
  }
  return run;
}

/// The space points of each track, innermost first.
std::vector<std::vector<std::uint32_t>> hitsPerTrack(
    const SpacePointContainer& spacePoints, const std::size_t nTracks) {
  auto trackColumn = spacePoints.column<std::uint32_t>("trackId");
  std::vector<std::vector<std::uint32_t>> hits(nTracks);
  for (const auto& sp : spacePoints) {
    hits.at(sp.extra(trackColumn)).push_back(sp.index());
  }
  for (auto& track : hits) {
    std::ranges::sort(track, {}, [&](const std::uint32_t index) {
      return spacePoints.at(index).r();
    });
  }
  return hits;
}

/// What the graph owed the tracks and what it delivered.
struct Completeness {
  std::size_t expectedEdges{};
  std::size_t foundEdges{};
  std::size_t expectedLinks{};
  std::size_t foundLinks{};
  std::string report;
};

/// Every consecutive triple of hits on a track is a grid triplet the search
/// owes it, so every consecutive pair is an edge of the graph and every
/// consecutive triple an edge connection between the two edges that share the
/// middle hit.
Completeness checkCompleteness(const GraphRun& run,
                               const SpacePointContainer& spacePoints,
                               const std::vector<Track>& tracks) {
  const std::vector<std::vector<std::uint32_t>> hits =
      hitsPerTrack(spacePoints, tracks.size());

  Completeness result;
  std::ostringstream missing;

  for (std::size_t track = 0; track < tracks.size(); ++track) {
    const std::vector<std::uint32_t>& own = hits[track];
    std::vector<std::optional<std::uint32_t>> edgeOf(own.size());

    for (std::size_t i = 0; i + 1 < own.size(); ++i) {
      ++result.expectedEdges;
      const std::pair<std::uint32_t, std::uint32_t> key{
          run.nodeOfSpacePoint[own[i]], run.nodeOfSpacePoint[own[i + 1]]};
      const auto found = run.edgeIndex.find(key);
      if (found == run.edgeIndex.end()) {
        missing << "  track " << track << " edge sp " << own[i] << "->"
                << own[i + 1] << " missing\n";
        continue;
      }
      ++result.foundEdges;
      edgeOf[i] = found->second;
    }

    for (std::size_t i = 0; i + 2 < own.size(); ++i) {
      ++result.expectedLinks;
      if (!edgeOf[i].has_value() || !edgeOf[i + 1].has_value()) {
        missing << "  track " << track << " link at sp " << own[i + 1]
                << " has no edge to link\n";
        continue;
      }
      // the outer edge records the inner one as a neighbour
      const std::uint32_t outer = *edgeOf[i + 1];
      const auto begin = run.neighbours.begin() + run.neighbourOffsets[outer];
      const auto end = run.neighbours.begin() + run.neighbourOffsets[outer + 1];
      if (std::find(begin, end, *edgeOf[i]) == end) {
        missing << "  track " << track << " link sp " << own[i] << "->"
                << own[i + 1] << "->" << own[i + 2] << " missing\n";
        continue;
      }
      ++result.foundLinks;
    }
  }

  result.report = missing.str();
  return result;
}

/// What the seeder returned, per track.
struct SeedTally {
  /// Seeds whose innermost hit belongs to this track.
  std::vector<std::size_t> seeds;
  /// Hits in the longest of them.
  std::vector<std::size_t> longest;
  /// Hits the track left in the container.
  std::vector<std::size_t> hits;
  /// Seeds mixing hits of more than one track.
  std::size_t impure{};
  /// The space point indices of every seed, in output order.
  std::vector<std::vector<std::uint32_t>> seedHits;
  std::string report;
};

/// Run the seeder and check the invariants that hold however dense the event
/// is: a seed holds hits of one track, innermost first.
SeedTally runSeeding(const SeederSetup& setup,
                     const SpacePointContainer& spacePoints,
                     const std::vector<Track>& tracks) {
  SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);
  const Experimental::TripletGbtsNodes nodes = makeNodes(setup, spacePoints);
  setup.selector.createSeeds(setup.seeder, nodes,
                             setup.graph->buildTheGraph(nodes, setup.bFieldInZ),
                             setup.filter, seeds);

  auto trackColumn = spacePoints.column<std::uint32_t>("trackId");

  SeedTally tally;
  tally.seeds.assign(tracks.size(), 0);
  tally.longest.assign(tracks.size(), 0);
  for (const std::vector<std::uint32_t>& own :
       hitsPerTrack(spacePoints, tracks.size())) {
    tally.hits.push_back(own.size());
  }

  for (const auto& seed : seeds) {
    const auto indices = seed.spacePointIndices();
    BOOST_REQUIRE_GE(indices.size(), 3u);
    tally.seedHits.emplace_back(indices.begin(), indices.end());
    const std::uint32_t track = spacePoints.at(indices[0]).extra(trackColumn);

    float previousR = -1.f;
    bool pure = true;
    for (const auto index : indices) {
      const auto sp = spacePoints.at(index);
      pure = pure && sp.extra(trackColumn) == track;
      BOOST_CHECK_GT(sp.r(), previousR);
      previousR = sp.r();
    }
    if (!pure) {
      ++tally.impure;
    }
    tally.seeds.at(track) += 1;
    tally.longest.at(track) = std::max(tally.longest.at(track), indices.size());
  }

  std::ostringstream detail;
  detail << "seeds: " << seeds.size() << " for " << tracks.size() << " tracks, "
         << tally.impure << " mixing tracks\n";
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    detail << "  track " << track << ": " << tally.seeds[track] << " seed(s)"
           << ", longest " << tally.longest[track] << " of "
           << tally.hits[track] << " hits\n";
  }
  tally.report = detail.str();
  return tally;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(TripletGbtsSeeding)

// Guards the fixture: the tracks have to actually be displaced, actually reach
// the layers, and actually stay inside the region the search is configured for.
BOOST_AUTO_TEST_CASE(DisplacedInputIsWellFormed) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDisplacedTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  float largestD0 = 0.f;
  float largestZ0 = 0.f;
  for (const Track& track : tracks) {
    largestD0 = std::max(largestD0, std::abs(track.d0));
    largestZ0 = std::max(largestZ0, std::abs(track.z0));
    BOOST_CHECK_LT(std::abs(track.d0), kD0Max);
    BOOST_CHECK_LT(std::abs(track.z0), kZ0Max);
  }
  // the region is covered, not just its middle
  BOOST_CHECK_GT(largestD0, 0.8f * kD0Max);
  BOOST_CHECK_GT(largestZ0, 0.8f * kZ0Max);

  // every track crosses every layer, so each owes two triplets
  BOOST_CHECK_EQUAL(spacePoints.size(), tracks.size() * detector.layers.size());
}

// The graph holds every edge and every edge connection the tracks owe it, out
// of the grid triplets alone, for tracks spread over the whole displaced region
// and for a tight bundle of them.
BOOST_AUTO_TEST_CASE(GraphFindsEveryEdgeAndConnection) {
  const ToyDetector detector = stripBarrelDetector();
  for (const bool dense : {false, true}) {
    const std::vector<Track> tracks =
        dense ? makeDenseDisplacedTracks() : makeDisplacedTracks();
    const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);
    const SeederSetup setup = makeSeeder(detector);
    const GraphRun run = buildGraph(setup, spacePoints);
    const Completeness result = checkCompleteness(run, spacePoints, tracks);

    BOOST_TEST_MESSAGE("dense " << dense << ": " << run.stats.triplets
                                << " triplets, " << run.nEdges << " edges, "
                                << run.nConnections << " connections\n"
                                << result.report);
    BOOST_CHECK_EQUAL(result.foundEdges, result.expectedEdges);
    BOOST_CHECK_EQUAL(result.foundLinks, result.expectedLinks);
    BOOST_CHECK_EQUAL(run.stats.tripletsOutsideTable, 0u);
  }
}

BOOST_AUTO_TEST_CASE(SeedsFromDisplacedTracks) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDisplacedTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const SeedTally tally = runSeeding(makeSeeder(detector), spacePoints, tracks);

  BOOST_TEST_MESSAGE(tally.report);

  // the tracks are far apart in phi, so clone removal should resolve the
  // branching to exactly one seed per track, holding every hit it left
  BOOST_CHECK_EQUAL(tally.impure, 0u);
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    BOOST_CHECK_EQUAL(tally.seeds[track], 1u);
    BOOST_CHECK_EQUAL(tally.longest[track], tally.hits[track]);
  }
}

// The same among a tight bundle of tracks, whose triplets the search mixes:
// every track still has to come back whole, which is the filter and clone
// removal resolving the branching.
BOOST_AUTO_TEST_CASE(SeedsFromDenseDisplacedTracks) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDenseDisplacedTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const SeedTally tally = runSeeding(makeSeeder(detector), spacePoints, tracks);

  BOOST_TEST_MESSAGE(tally.report);

  BOOST_CHECK_EQUAL(tally.impure, 0u);
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    BOOST_CHECK_GE(tally.seeds[track], 1u);
    BOOST_CHECK_EQUAL(tally.longest[track], tally.hits[track]);
  }
}

// Pairwise clone removal calls a candidate a clone only where the default
// would too, so it keeps every seed the default keeps.
BOOST_AUTO_TEST_CASE(PairwiseCloneRemovalKeepsEverySeedTheDefaultKeeps) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDenseDisplacedTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const SeedTally byDefault =
      runSeeding(makeSeeder(detector), spacePoints, tracks);

  for (const bool dropFullyTaken : {false, true}) {
    const SeederSetup setup = makeSeeder(
        detector, {},
        [dropFullyTaken](Experimental::GraphBasedTrackSeeder::Config& /*cfg*/,
                         Experimental::TripletGbtsSeedSelector::Config& sel) {
          sel.pairwiseCloneRemoval = true;
          sel.dropFullyTakenTriplets = dropFullyTaken;
        });
    const SeedTally pairwise = runSeeding(setup, spacePoints, tracks);

    BOOST_TEST_MESSAGE("dropFullyTakenTriplets " << dropFullyTaken << ": "
                                                 << pairwise.report);

    std::multiset<std::vector<std::uint32_t>> kept(pairwise.seedHits.begin(),
                                                   pairwise.seedHits.end());
    for (const std::vector<std::uint32_t>& seed : byDefault.seedHits) {
      const auto found = kept.find(seed);
      BOOST_CHECK(found != kept.end());
      if (found != kept.end()) {
        kept.erase(found);
      }
    }
  }
}

// Dropping the triplets that share a hit with a longer kept seed removes seeds
// and still leaves every track whole.
BOOST_AUTO_TEST_CASE(DropTripletsSharingLongerKeepsEveryTrackWhole) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDenseDisplacedTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const auto pairwise = [&](const bool dropSharingLonger) {
    return runSeeding(
        makeSeeder(detector, {},
                   [dropSharingLonger](
                       Experimental::GraphBasedTrackSeeder::Config& /*cfg*/,
                       Experimental::TripletGbtsSeedSelector::Config& sel) {
                     sel.pairwiseCloneRemoval = true;
                     sel.dropTripletsSharingLonger = dropSharingLonger;
                   }),
        spacePoints, tracks);
  };
  const SeedTally alone = pairwise(false);
  const SeedTally dropping = pairwise(true);

  BOOST_TEST_MESSAGE("pairwise alone: " << alone.report);
  BOOST_TEST_MESSAGE("dropping triplets sharing longer: " << dropping.report);

  BOOST_CHECK_LE(dropping.seedHits.size(), alone.seedHits.size());
  BOOST_CHECK_LE(dropping.impure, alone.impure);
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    BOOST_CHECK_GE(dropping.seeds[track], 1u);
    BOOST_CHECK_EQUAL(dropping.longest[track], dropping.hits[track]);
  }
}

// Capping the three space point candidates per middle node leaves at most that
// many triplets on each middle space point, and still every track whole.
BOOST_AUTO_TEST_CASE(TripletCapPerMiddleNodeKeepsEveryTrackWhole) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDenseDisplacedTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  for (const std::uint32_t cap : {1u, 2u}) {
    const SeedTally tally = runSeeding(
        makeSeeder(detector, {},
                   [cap](Experimental::GraphBasedTrackSeeder::Config& cfg,
                         Experimental::TripletGbtsSeedSelector::Config& sel) {
                     sel.pairwiseCloneRemoval = true;
                     sel.maxTripletsPerMiddleNode = cap;
                     // a split seed's halves are triplets the cap never saw
                     cfg.maxSeedSplitEta = 0.f;
                   }),
        spacePoints, tracks);

    BOOST_TEST_MESSAGE("cap " << cap << ": " << tally.report);

    std::map<std::uint32_t, std::uint32_t> perMiddle;
    for (const std::vector<std::uint32_t>& seed : tally.seedHits) {
      if (seed.size() == 3) {
        ++perMiddle[seed[1]];
      }
    }
    for (const auto& [middle, n] : perMiddle) {
      BOOST_CHECK_LE(n, cap);
    }
    for (std::size_t track = 0; track < tracks.size(); ++track) {
      BOOST_CHECK_GE(tally.seeds[track], 1u);
      BOOST_CHECK_EQUAL(tally.longest[track], tally.hits[track]);
    }
  }
}

// With stripLayersFromTable the layer connection table decides which layers
// the graph links: every edge then joins a pair of layers the table connects,
// and taking a connection out of the table takes its edges out of the graph.
BOOST_AUTO_TEST_CASE(StripLayersFromTableKeepToTheTable) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDisplacedTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);
  auto layerColumn =
      spacePoints.column<Experimental::GbtsLayerIndex>("gbtsLayerIndex");

  // the dense layer index of each table connection, outer layer first
  std::set<std::pair<std::uint32_t, std::uint32_t>> table;
  for (std::size_t i = 0; i + 1 < detector.layers.size(); ++i) {
    table.emplace(static_cast<std::uint32_t>(i + 1),
                  static_cast<std::uint32_t>(i));
  }
  const auto layerPairs = [&](const GraphRun& run) {
    std::set<std::pair<std::uint32_t, std::uint32_t>> pairs;
    for (const auto& edge : run.edges) {
      pairs.emplace(
          spacePoints.at(run.spacePointOfNode[edge.n2]).extra(layerColumn),
          spacePoints.at(run.spacePointOfNode[edge.n1]).extra(layerColumn));
    }
    return pairs;
  };
  const auto withTable = [](Experimental::TripletGbtsGraphBuilder::Config& c) {
    c.stripLayersFromTable = true;
  };

  // the full table: every edge on a connection, every track still whole
  const SeederSetup full = makeSeeder(detector, withTable);
  const GraphRun fullRun = buildGraph(full, spacePoints);
  BOOST_TEST_MESSAGE("full table: " << fullRun.stats.triplets << " triplets, "
                                    << fullRun.stats.tripletsOutsideTable
                                    << " outside it");
  for (const auto& pair : layerPairs(fullRun)) {
    BOOST_CHECK(table.contains(pair));
  }
  const Completeness complete =
      checkCompleteness(fullRun, spacePoints, tracks);
  BOOST_CHECK_EQUAL(complete.foundLinks, complete.expectedLinks);

  // the middle connection taken out: none of its edges left
  const SeederSetup cut = makeSeeder(
      detector, withTable, {},
      [](std::vector<Experimental::GbtsLayerConnection>& connections) {
        std::erase_if(connections,
                      [](const Experimental::GbtsLayerConnection& c) {
                        return c.src == 82000 && c.dst == 81000;
                      });
      });
  const GraphRun cutRun = buildGraph(cut, spacePoints);
  BOOST_CHECK_EQUAL(cutRun.stats.triplets, fullRun.stats.triplets);
  BOOST_CHECK_GT(cutRun.stats.tripletsOutsideTable,
                 fullRun.stats.tripletsOutsideTable);
  BOOST_CHECK(!layerPairs(cutRun).contains({2u, 1u}));
}

// With stripSearchFromTable the grid search itself combines only the layers
// the table connects: no triplet is left for stripLayersFromTable to drop,
// every edge joins a connected pair, the tracks stay whole, and a connection
// taken out of the table is never searched.
BOOST_AUTO_TEST_CASE(StripSearchFromTableSearchesOnlyTheTable) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDisplacedTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);
  auto layerColumn =
      spacePoints.column<Experimental::GbtsLayerIndex>("gbtsLayerIndex");

  // the dense layer index of each table connection, outer layer first
  std::set<std::pair<std::uint32_t, std::uint32_t>> table;
  for (std::size_t i = 0; i + 1 < detector.layers.size(); ++i) {
    table.emplace(static_cast<std::uint32_t>(i + 1),
                  static_cast<std::uint32_t>(i));
  }
  const auto layerPairs = [&](const GraphRun& run) {
    std::set<std::pair<std::uint32_t, std::uint32_t>> pairs;
    for (const auto& edge : run.edges) {
      pairs.emplace(
          spacePoints.at(run.spacePointOfNode[edge.n2]).extra(layerColumn),
          spacePoints.at(run.spacePointOfNode[edge.n1]).extra(layerColumn));
    }
    return pairs;
  };
  const auto withTable = [](Experimental::TripletGbtsGraphBuilder::Config& c) {
    c.stripLayersFromTable = true;
  };
  const auto withSearch = [](Experimental::TripletGbtsGraphBuilder::Config& c) {
    c.stripLayersFromTable = true;
    c.stripSearchFromTable = true;
  };

  const GraphRun filteredRun =
      buildGraph(makeSeeder(detector, withTable), spacePoints);
  const SeederSetup search = makeSeeder(detector, withSearch);
  const GraphRun searchRun = buildGraph(search, spacePoints);
  BOOST_TEST_MESSAGE("filtered: " << filteredRun.stats.triplets
                                  << " triplets, "
                                  << filteredRun.stats.tripletsOutsideTable
                                  << " outside the table; searched: "
                                  << searchRun.stats.triplets << " triplets");
  BOOST_CHECK_EQUAL(searchRun.stats.tripletsOutsideTable, 0u);
  BOOST_CHECK_LE(searchRun.stats.triplets, filteredRun.stats.triplets);
  for (const auto& pair : layerPairs(searchRun)) {
    BOOST_CHECK(table.contains(pair));
  }
  const Completeness complete =
      checkCompleteness(searchRun, spacePoints, tracks);
  BOOST_CHECK_EQUAL(complete.foundLinks, complete.expectedLinks);
  BOOST_TEST_MESSAGE(complete.report);

  const SeedTally tally = runSeeding(search, spacePoints, tracks);
  BOOST_TEST_MESSAGE(tally.report);
  BOOST_CHECK_EQUAL(tally.impure, 0u);
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    BOOST_CHECK_GE(tally.seeds[track], 1u);
    BOOST_CHECK_EQUAL(tally.longest[track], tally.hits[track]);
  }

  // the middle connection taken out: no triplet uses it, so fewer are made
  const SeederSetup cut = makeSeeder(
      detector, withSearch, {},
      [](std::vector<Experimental::GbtsLayerConnection>& connections) {
        std::erase_if(connections,
                      [](const Experimental::GbtsLayerConnection& c) {
                        return c.src == 82000 && c.dst == 81000;
                      });
      });
  const GraphRun cutRun = buildGraph(cut, spacePoints);
  BOOST_CHECK_EQUAL(cutRun.stats.tripletsOutsideTable, 0u);
  BOOST_CHECK_LT(cutRun.stats.triplets, searchRun.stats.triplets);
  BOOST_CHECK(!layerPairs(cutRun).contains({2u, 1u}));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
