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
#include "Acts/Seeding/DisplacedGbtsGraph.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsLayerConnection.hpp"
#include "Acts/Seeding/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"
#include "Acts/Seeding/detail/DisplacedGraphTypes.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Efficiency harness for the displaced (large radius) GBTS graph. Unlike the
// prompt suite next door this pins nothing recorded from the implementation:
// it generates displaced helices, works out from the geometry alone which
// edges and which edge connections the graph owes them, and checks the graph
// produced every one.

namespace Acts::Test {

namespace {

using namespace Acts::UnitLiterals;
using Experimental::GbtsLayerTechnology;
using Experimental::GbtsLayerType;
using Experimental::detail::DisplacedGbtsEdge;

/// The region the displaced graph is asked to cover.
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

std::shared_ptr<Experimental::GbtsGeometry> makeGeometry(
    const ToyDetector& detector) {
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
    layers.push_back(layer);
  }

  std::vector<Experimental::GbtsLayerConnection> connections;
  connections.reserve(detector.links.size());
  for (const auto& [src, dst] : detector.links) {
    connections.push_back({src, dst});
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
/// the one thing the displaced graph may not assume, so these carry a real
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
    const float a = static_cast<float>(sense) * t + 0.5f * std::numbers::pi_v<float>;
    const float dx = static_cast<float>(sense) * (vx * std::cos(a) - vy * std::sin(a));
    const float dy = static_cast<float>(sense) * (vx * std::sin(a) + vy * std::cos(a));
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
      SpacePointColumns::Phi | SpacePointColumns::StripCalibrationDetails);

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

/// The seeder and everything the displaced graph needs.
struct SeederSetup {
  Experimental::GraphBasedTrackSeeder seeder;
  Experimental::DisplacedGbtsGraph graph;
  Experimental::GbtsTrackingFilter filter;
  Experimental::GbtsRoiDescriptor roi;
  Experimental::GraphBasedTrackSeeder::Options options;
};

SeederSetup makeSeeder(const ToyDetector& detector,
                       const bool calibrateStrips = true,
                       const bool matchBeforeCreate = false) {
  auto geometry = makeGeometry(detector);

  const auto makeLogger = []() -> std::unique_ptr<const Logger> {
    return getDefaultLogger("DisplacedGbtsTest", Logging::Level::WARNING);
  };

  Experimental::GraphBasedTrackSeeder::Config config;
  // the toy setup has no tau lookup table and no cluster widths
  config.useClusterWidthCuts = false;

  Experimental::GbtsGraphConfig graphConfig;
  graphConfig.minPt = 1_GeV;
  graphConfig.d0Max = kD0Max;
  graphConfig.minZ0 = -kZ0Max;
  graphConfig.maxZ0 = kZ0Max;
  graphConfig.maxOuterRadius = kMaxOuterRadius;
  graphConfig.calibrateStrips = calibrateStrips;
  graphConfig.matchBeforeCreate = matchBeforeCreate;
  // a link here already cost a triplet fit, so three nodes are a seed
  graphConfig.minSeedLevel = 2;

  // The filter walks the chains the graph found and has its own beamspot
  // assumption, defaulting to 170mm. Left there it would throw away the
  // displaced chains the graph was configured to go and find.
  Experimental::GbtsTrackingFilter::Config filterConfig;
  filterConfig.maxZ0 = kZ0Max;

  return SeederSetup{
      .seeder = Experimental::GraphBasedTrackSeeder(
          Experimental::GraphBasedTrackSeeder::DerivedConfig(config), geometry,
          makeLogger()),
      .graph = Experimental::DisplacedGbtsGraph(graphConfig, geometry,
                                                makeLogger()),
      .filter = Experimental::GbtsTrackingFilter(filterConfig, geometry,
                                                makeLogger()),
      .roi = Experimental::GbtsRoiDescriptor(-4.5, 4.5, -kZ0Max, kZ0Max),
      .options = Experimental::GraphBasedTrackSeeder::Options{.bFieldInZ = 2_T},
  };
}

/// The graph, plus what is needed to read it back against the tracks.
struct GraphRun {
  std::vector<DisplacedGbtsEdge> edges;
  std::uint32_t nEdges{};
  std::uint32_t nConnections{};
  /// Node index of each space point, by space point index.
  std::vector<std::uint32_t> nodeOfSpacePoint;
  /// Edge index by (inner node, outer node).
  std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> edgeIndex;
};

GraphRun buildGraph(const SeederSetup& setup,
                    const SpacePointContainer& spacePoints) {
  Experimental::GbtsNodeStorage storage = setup.seeder.makeNodeStorage();

  auto layerColumn =
      spacePoints.column<Experimental::GbtsLayerIndex>("gbtsLayerIndex");
  auto clusterWidthColumn = spacePoints.column<float>("clusterWidth");
  auto localPositionColumn = spacePoints.column<float>("localPositionY");
  storage.extend(spacePoints, layerColumn, clusterWidthColumn,
                 localPositionColumn);
  storage.finalize();

  GraphRun run;
  const auto stats = setup.graph.buildTheGraph(
      setup.roi, storage, run.edges, setup.options.bFieldInZ);
  run.nEdges = stats.first;
  run.nConnections = stats.second;

  run.nodeOfSpacePoint.assign(spacePoints.size(), 0);
  for (std::uint32_t node = 0; node < storage.numberOfNodes(); ++node) {
    run.nodeOfSpacePoint.at(storage.spacePointIndex(node)) = node;
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

/// Every consecutive pair of hits on a track is an edge the graph owes it, and
/// every consecutive triple an edge connection between the two edges that
/// share the middle hit. The layers are linked only to their neighbours, so
/// this is the whole of what the graph should find.
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
      const DisplacedGbtsEdge& outer = run.edges[*edgeOf[i + 1]];
      const auto begin = outer.vNei.begin();
      const auto end = begin + outer.nNei;
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
  std::string report;
};

/// Run the seeder and check the invariants that hold however dense the event
/// is: a seed holds hits of one track, innermost first.
SeedTally runSeeding(const SeederSetup& setup,
                     const SpacePointContainer& spacePoints,
                     const std::vector<Track>& tracks) {
  SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);
  setup.seeder.createSeeds(spacePoints, setup.roi, setup.graph, setup.filter,
                           setup.options, seeds);

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
  detail << "seeds: " << seeds.size() << " for " << tracks.size()
         << " tracks, " << tally.impure << " mixing tracks\n";
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    detail << "  track " << track << ": " << tally.seeds[track] << " seed(s)"
           << ", longest " << tally.longest[track] << " of "
           << tally.hits[track] << " hits\n";
  }
  tally.report = detail.str();
  return tally;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(DisplacedGbtsSeeding)

// Guards the fixture: the tracks have to actually be displaced, actually reach
// the layers, and actually stay inside the region the graph is configured for.
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

  // every track crosses every layer, so each owes three edges and two links
  BOOST_CHECK_EQUAL(spacePoints.size(), tracks.size() * detector.layers.size());

  auto layerColumn =
      spacePoints.column<Experimental::GbtsLayerIndex>("gbtsLayerIndex");
  for (const auto& sp : spacePoints) {
    BOOST_CHECK_LT(sp.extra(layerColumn), detector.layers.size());
    BOOST_CHECK_LE(std::abs(sp.z()), kBarrelHalfZ);
    BOOST_CHECK_GT(sp.r(), kD0Max);
    BOOST_CHECK_LT(sp.r(), kMaxOuterRadius);
  }
}

// The point of the suite: on hits that came from real displaced helices the
// graph must find every doublet and every doublet pair they imply.
BOOST_AUTO_TEST_CASE(GraphFindsEveryEdgeAndConnection) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDisplacedTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const SeederSetup setup = makeSeeder(detector);
  const GraphRun run = buildGraph(setup, spacePoints);

  const Completeness result = checkCompleteness(run, spacePoints, tracks);

  BOOST_TEST_MESSAGE("graph: " << run.nEdges << " edges, " << run.nConnections
                               << " connections");
  BOOST_TEST_MESSAGE("owed: " << result.expectedEdges << " edges, "
                              << result.expectedLinks << " connections");
  if (!result.report.empty()) {
    BOOST_TEST_MESSAGE("missing:\n" << result.report);
  }

  BOOST_CHECK_EQUAL(result.foundEdges, result.expectedEdges);
  BOOST_CHECK_EQUAL(result.foundLinks, result.expectedLinks);
}

// The same, with the strip points walked off along their own strips, so the
// graph only finds them by resolving each triplet against its own fit.
BOOST_AUTO_TEST_CASE(GraphFindsEveryEdgeAfterResolvingTheStrips) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDisplacedTracks();
  constexpr float kWalk = 8.f;
  const SpacePointContainer spacePoints =
      makeSpacePoints(detector, tracks, kWalk);

  const SeederSetup setup = makeSeeder(detector, /*calibrateStrips=*/true);
  const GraphRun run = buildGraph(setup, spacePoints);

  const Completeness result = checkCompleteness(run, spacePoints, tracks);

  BOOST_TEST_MESSAGE("walked graph: " << run.nEdges << " edges, "
                                      << run.nConnections << " connections");
  if (!result.report.empty()) {
    BOOST_TEST_MESSAGE("missing:\n" << result.report);
  }

  BOOST_CHECK_EQUAL(result.foundEdges, result.expectedEdges);
  BOOST_CHECK_EQUAL(result.foundLinks, result.expectedLinks);
}

// With the strip resolution switched off nothing moves, so the geometric
// precut on the turn between two edges is allowed to run. It is the one cut
// that can drop a candidate before the fit ever sees it, so the tracks have to
// survive it whole.
BOOST_AUTO_TEST_CASE(GraphFindsEveryEdgeWithTheTurnPrecutActive) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDisplacedTracks();
  // unwalked, so the uncalibrated points are already the true crossings
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const SeederSetup setup = makeSeeder(detector, /*calibrateStrips=*/false);
  const GraphRun run = buildGraph(setup, spacePoints);

  const Completeness result = checkCompleteness(run, spacePoints, tracks);

  BOOST_TEST_MESSAGE("uncalibrated graph: " << run.nEdges << " edges, "
                                            << run.nConnections
                                            << " connections");
  if (!result.report.empty()) {
    BOOST_TEST_MESSAGE("missing:\n" << result.report);
  }

  BOOST_CHECK_EQUAL(result.foundEdges, result.expectedEdges);
  BOOST_CHECK_EQUAL(result.foundLinks, result.expectedLinks);
}

// Every barrel layer is given a depth counting outwards from the innermost,
// whatever it is made of. `barrelOrder` cannot serve: it numbers the pixel
// barrel alone, so on this detector it is -1 throughout and anything keyed on
// it never fires.
BOOST_AUTO_TEST_CASE(BarrelLayersAreGivenADepth) {
  const ToyDetector detector = stripBarrelDetector();
  const auto geometry = makeGeometry(detector);

  BOOST_REQUIRE_EQUAL(geometry->numLayers(), detector.layers.size());

  for (Experimental::GbtsLayerIndex i = 0; i < geometry->numLayers(); ++i) {
    const Experimental::GbtsLayerDescription& layer =
        geometry->layerDescription(i);

    // the toy layers are handed over innermost first, so depth follows
    BOOST_CHECK_EQUAL(layer.depth, static_cast<std::int32_t>(i));
    // and the thing it replaces is blind to them
    BOOST_CHECK_EQUAL(layer.barrelOrder, -1);
  }
}

// With a depth to key on, matchBeforeCreate reaches this detector at all. It
// is a filter on edge creation, so the tracks have to survive it whole.
//
// The dense fixture, since it only filters a node carrying more than
// `matchBeforeCreateMaxEdges` outgoing edges and well separated tracks never
// give one that many.
BOOST_AUTO_TEST_CASE(MatchBeforeCreateKeepsEveryEdge) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDenseDisplacedTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const SeederSetup off = makeSeeder(detector, true, /*matchBeforeCreate=*/false);
  const SeederSetup on = makeSeeder(detector, true, /*matchBeforeCreate=*/true);

  const GraphRun without = buildGraph(off, spacePoints);
  const GraphRun with = buildGraph(on, spacePoints);

  BOOST_TEST_MESSAGE("matchBeforeCreate off: " << without.nEdges << " edges, "
                                               << without.nConnections
                                               << " connections");
  BOOST_TEST_MESSAGE("matchBeforeCreate on : " << with.nEdges << " edges, "
                                               << with.nConnections
                                               << " connections");

  // it has to actually engage, or the test below proves nothing
  BOOST_CHECK_LT(with.nEdges, without.nEdges);

  // and every edge and connection the tracks owe is still there
  const Completeness result = checkCompleteness(with, spacePoints, tracks);
  if (!result.report.empty()) {
    BOOST_TEST_MESSAGE("missing:\n" << result.report);
  }
  BOOST_CHECK_EQUAL(result.foundEdges, result.expectedEdges);
  BOOST_CHECK_EQUAL(result.foundLinks, result.expectedLinks);
}

// Guards the dense fixture: it is only worth running if the tracks really do
// overlap, which is measured by how much work the graph is made to do beyond
// the edges the tracks themselves owe it.
BOOST_AUTO_TEST_CASE(DenseInputIsActuallyDense) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDenseDisplacedTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  BOOST_CHECK_EQUAL(spacePoints.size(), tracks.size() * detector.layers.size());

  // the hits of one layer, spread over a hand's width of azimuth at most
  auto layerColumn =
      spacePoints.column<Experimental::GbtsLayerIndex>("gbtsLayerIndex");
  float minPhi = std::numbers::pi_v<float>;
  float maxPhi = -std::numbers::pi_v<float>;
  for (const auto& sp : spacePoints) {
    if (sp.extra(layerColumn) != 0) {
      continue;
    }
    minPhi = std::min(minPhi, sp.phi());
    maxPhi = std::max(maxPhi, sp.phi());
  }
  BOOST_TEST_MESSAGE("innermost layer spans " << (maxPhi - minPhi)
                                              << " rad over " << tracks.size()
                                              << " tracks");
  BOOST_CHECK_LT(maxPhi - minPhi, 0.2f);
}

// The same completeness demand under real combinatorics: the tracks overlap,
// so the graph builds far more edges than they owe it and has to keep all of
// theirs among them.
BOOST_AUTO_TEST_CASE(GraphFindsEveryEdgeAmongDenseTracks) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDenseDisplacedTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const SeederSetup setup = makeSeeder(detector);
  const GraphRun run = buildGraph(setup, spacePoints);

  const Completeness result = checkCompleteness(run, spacePoints, tracks);

  BOOST_TEST_MESSAGE("dense graph: "
                     << run.nEdges << " edges (" << result.expectedEdges
                     << " owed), " << run.nConnections << " connections ("
                     << result.expectedLinks << " owed)");
  if (!result.report.empty()) {
    BOOST_TEST_MESSAGE("missing:\n" << result.report);
  }

  BOOST_CHECK_EQUAL(result.foundEdges, result.expectedEdges);
  BOOST_CHECK_EQUAL(result.foundLinks, result.expectedLinks);

  // the point of the fixture: the graph is carrying real combinatorics here,
  // not just the tracks' own edges
  BOOST_CHECK_GT(run.nEdges, 2 * result.expectedEdges);
}

// Everything above feeds seeds: one per track, holding all four hits.
BOOST_AUTO_TEST_CASE(SeedsFromDisplacedTracks) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDisplacedTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const SeederSetup setup = makeSeeder(detector);
  const SeedTally tally = runSeeding(setup, spacePoints, tracks);

  BOOST_TEST_MESSAGE(tally.report);

  // the tracks are far apart in phi, so clone removal should resolve the
  // branching to exactly one seed per track, holding every hit it left
  BOOST_CHECK_EQUAL(tally.impure, 0u);
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    BOOST_CHECK_EQUAL(tally.seeds[track], 1u);
    BOOST_CHECK_EQUAL(tally.longest[track], tally.hits[track]);
  }
}

// The same, out of a graph carrying thirty times the edges the tracks owe it.
// Every track still has to come back whole, which is the filter and clone
// removal resolving the branching rather than the graph handing them an easy
// problem.
BOOST_AUTO_TEST_CASE(SeedsFromDenseDisplacedTracks) {
  const ToyDetector detector = stripBarrelDetector();
  const std::vector<Track> tracks = makeDenseDisplacedTracks();
  const SpacePointContainer spacePoints = makeSpacePoints(detector, tracks);

  const SeederSetup setup = makeSeeder(detector);
  const SeedTally tally = runSeeding(setup, spacePoints, tracks);

  BOOST_TEST_MESSAGE(tally.report);

  // every track back whole, and no seed built out of two of them
  BOOST_CHECK_EQUAL(tally.impure, 0u);
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    BOOST_CHECK_GE(tally.seeds[track], 1u);
    BOOST_CHECK_EQUAL(tally.longest[track], tally.hits[track]);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
