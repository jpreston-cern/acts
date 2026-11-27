// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/TripletGbtsGraphBuilder.hpp"

#include "Acts/Seeding/GbtsGeometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <functional>
#include <numbers>
#include <queue>
#include <span>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Acts::Experimental {

namespace {

/// Signed transverse curvature of the circle through three nodes (1/mm)
float curvature(const std::array<float, 4>& a, const std::array<float, 4>& b,
                const std::array<float, 4>& c) {
  const float cross =
      (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
  const float d12 = std::hypot(b[0] - a[0], b[1] - a[1]);
  const float d23 = std::hypot(c[0] - b[0], c[1] - b[1]);
  const float d13 = std::hypot(c[0] - a[0], c[1] - a[1]);
  const float denominator = d12 * d23 * d13;
  return denominator > 0.f ? 2.f * cross / denominator : 0.f;
}

/// dz/dr from the inner to the outer node
float tau(const std::array<float, 4>& inner,
          const std::array<float, 4>& outer) {
  const float dr = outer[3] - inner[3];
  return dr != 0.f ? (outer[2] - inner[2]) / dr : 0.f;
}

}  // namespace

template <typename edge_of_t>
void TripletGbtsGraphBuilder::extendIntoPixels(
    const detail::GbtsNodeView& nodeView, const std::uint32_t nNodes,
    std::vector<EdgeType>& edges, edge_of_t& edgeOf,
    std::vector<std::tuple<std::uint32_t, float, std::uint32_t>>& rawLinks,
    std::span<const std::uint8_t> allowedPairs, const std::uint32_t nLayers,
    Stats& stats) const {
  const auto& positions = nodeView.positions;
  const auto isPixel = [&](const std::uint32_t n) {
    return nodeView.strip(n) == nullptr;
  };

  // the pixel nodes per layer, in azimuth bins, each sorted by z on a barrel
  // layer and by r on an endcap layer, so that an edge visits only the nodes
  // inside a window its dz/dr tolerance allows; the exact cuts follow
  constexpr std::uint32_t kPhiBins = 64;
  const float pi = std::numbers::pi_v<float>;
  const auto phiBin = [&](const float phi) {
    const auto b =
        static_cast<std::int32_t>((phi + pi) / (2.f * pi) * kPhiBins);
    return static_cast<std::uint32_t>(
        std::clamp<std::int32_t>(b, 0, kPhiBins - 1));
  };
  struct PixelLayer {
    float rMin = std::numeric_limits<float>::max();
    float rMax = std::numeric_limits<float>::lowest();
    float zMin = std::numeric_limits<float>::max();
    float zMax = std::numeric_limits<float>::lowest();
    bool barrel = true;
    // the layer's dense index
    std::uint32_t index = 0;
    // per azimuth bin: (z or r, azimuth, node)
    std::array<std::vector<std::tuple<float, float, std::uint32_t>>, kPhiBins>
        bins;
  };
  std::vector<PixelLayer> layers;
  std::unordered_map<std::uint32_t, std::size_t> layerOf;
  for (std::uint32_t n = 0; n < nNodes; ++n) {
    if (!isPixel(n)) {
      continue;
    }
    const auto [it, inserted] = layerOf.try_emplace(
        static_cast<std::uint32_t>(nodeView.layers[n]), layers.size());
    if (inserted) {
      layers.emplace_back();
      layers.back().index = static_cast<std::uint32_t>(nodeView.layers[n]);
    }
    PixelLayer& layer = layers[it->second];
    layer.rMin = std::min(layer.rMin, positions[n][3]);
    layer.rMax = std::max(layer.rMax, positions[n][3]);
    layer.zMin = std::min(layer.zMin, positions[n][2]);
    layer.zMax = std::max(layer.zMax, positions[n][2]);
  }
  if (layers.empty()) {
    return;
  }
  for (PixelLayer& layer : layers) {
    layer.barrel = (layer.rMax - layer.rMin) < (layer.zMax - layer.zMin);
  }
  for (std::uint32_t n = 0; n < nNodes; ++n) {
    if (!isPixel(n)) {
      continue;
    }
    PixelLayer& layer =
        layers[layerOf.at(static_cast<std::uint32_t>(nodeView.layers[n]))];
    const float phi = std::atan2(positions[n][1], positions[n][0]);
    layer.bins[phiBin(phi)].emplace_back(
        layer.barrel ? positions[n][2] : positions[n][3], phi, n);
  }
  for (PixelLayer& layer : layers) {
    for (auto& bin : layer.bins) {
      std::ranges::sort(bin);
    }
  }

  // with pixelLinksFromTable a strip layer is extended when the table connects
  // it to a pixel layer
  const bool linksFromTable =
      m_cfg.pixelLinksFromTable && !allowedPairs.empty();
  const float minDeltaR = m_cfg.pixelMinDeltaR;
  std::vector<std::uint8_t> reachesPixels;
  if (linksFromTable) {
    reachesPixels.assign(nLayers, 0);
    for (std::uint32_t outer = 0; outer < nLayers; ++outer) {
      for (const PixelLayer& layer : layers) {
        if (allowedPairs[static_cast<std::size_t>(outer) * nLayers +
                         layer.index] != 0) {
          reachesPixels[outer] = 1;
          break;
        }
      }
    }
  }

  // the curvature range of the triplets each edge is the inner edge of
  constexpr float kNone = std::numeric_limits<float>::quiet_NaN();
  std::vector<std::pair<float, float>> curvatureRange(edges.size(),
                                                      {kNone, kNone});
  const auto widen = [&](const std::uint32_t e, const float k) {
    if (e >= curvatureRange.size()) {
      curvatureRange.resize(e + 1, {kNone, kNone});
    }
    auto& [lo, hi] = curvatureRange[e];
    if (std::isnan(lo)) {
      lo = hi = k;
    } else {
      lo = std::min(lo, k);
      hi = std::max(hi, k);
    }
  };
  for (const auto& [outer, negQuality, inner] : rawLinks) {
    widen(inner,
          curvature(positions[edges[inner].n1], positions[edges[inner].n2],
                    positions[edges[outer].n2]));
  }

  // the pixel nodes within the azimuth window of node `a` and, per layer the
  // connection table allows from the layer of `a` (if it is used), inside the
  // z (barrel) or r (endcap) range that a dz/dr from the node in
  // [tauLo, tauHi] and a radius below r(a) - pixelMinDeltaR allow
  const auto forPixelsNear = [&](const std::uint32_t a, const float tauLo,
                                 const float tauHi, const float w,
                                 auto&& visit) {
    const float ra = positions[a][3];
    const float za = positions[a][2];
    const float phi = std::atan2(positions[a][1], positions[a][0]);
    const float rLimit = ra - minDeltaR;
    constexpr float kMargin = 1e-2f;
    const std::uint32_t layerA = nodeView.layers[a];
    for (const PixelLayer& layer : layers) {
      if (layer.rMin > rLimit) {
        continue;
      }
      if (!allowedPairs.empty() &&
          allowedPairs[static_cast<std::size_t>(layerA) * nLayers +
                       layer.index] == 0) {
        continue;
      }
      float keyLo = std::numeric_limits<float>::lowest();
      float keyHi = std::numeric_limits<float>::max();
      if (layer.barrel) {
        // z = za - tau * (ra - r), bilinear in tau and r: extremes at corners
        const float rHi = std::min(layer.rMax, rLimit);
        keyLo = std::numeric_limits<float>::max();
        keyHi = std::numeric_limits<float>::lowest();
        for (const float t : {tauLo, tauHi}) {
          for (const float r : {layer.rMin, rHi}) {
            const float z = za - t * (ra - r);
            keyLo = std::min(keyLo, z);
            keyHi = std::max(keyHi, z);
          }
        }
      } else if (tauLo > 0.f || tauHi < 0.f) {
        // r = ra - (za - z) / tau, monotonic in z and in tau of one sign
        keyLo = std::numeric_limits<float>::max();
        keyHi = std::numeric_limits<float>::lowest();
        for (const float t : {tauLo, tauHi}) {
          for (const float z : {layer.zMin, layer.zMax}) {
            const float r = ra - (za - z) / t;
            keyLo = std::min(keyLo, r);
            keyHi = std::max(keyHi, r);
          }
        }
        keyHi = std::min(keyHi, rLimit);
      } else {
        // a dz/dr window around 0: |za - z| <= max|tau| * (ra - r) bounds r by
        // ra less the distance in z from the node to the layer over max|tau|
        const float dz = za < layer.zMin   ? layer.zMin - za
                         : za > layer.zMax ? za - layer.zMax
                                           : 0.f;
        const float maxAbsTau = std::max(std::abs(tauLo), std::abs(tauHi));
        if (maxAbsTau > 0.f) {
          keyHi = std::min(rLimit, ra - dz / maxAbsTau);
        } else {
          keyHi = dz > 0.f ? std::numeric_limits<float>::lowest() : rLimit;
        }
      }
      keyLo -= kMargin;
      keyHi += kMargin;
      // a layer the window misses altogether
      if (layer.barrel ? (keyHi < layer.zMin || keyLo > layer.zMax)
                       : (keyHi < layer.rMin || keyLo > layer.rMax)) {
        continue;
      }
      const auto scanBins = [&](const std::uint32_t firstBin,
                                const std::uint32_t lastBin) {
        for (std::uint32_t i = firstBin; i <= lastBin; ++i) {
          const auto& bin = layer.bins[i];
          auto it = std::ranges::lower_bound(
              bin, keyLo, {}, [](const auto& e) { return std::get<0>(e); });
          for (; it != bin.end() && std::get<0>(*it) <= keyHi; ++it) {
            float dPhi = std::get<1>(*it) - phi;
            if (dPhi > pi) {
              dPhi -= 2.f * pi;
            } else if (dPhi < -pi) {
              dPhi += 2.f * pi;
            }
            if (std::abs(dPhi) <= w) {
              visit(std::get<2>(*it));
            }
          }
        }
      };
      // only the azimuth bins the window overlaps (the binning is monotonic
      // in azimuth, so a node inside the window is in one of them)
      const float lo = phi - w;
      const float hi = phi + w;
      if (2.f * w >= 2.f * pi) {
        scanBins(0, kPhiBins - 1);
      } else if (lo < -pi) {
        scanBins(phiBin(lo + 2.f * pi), kPhiBins - 1);
        scanBins(0, phiBin(hi));
      } else if (hi > pi) {
        scanBins(phiBin(lo), kPhiBins - 1);
        scanBins(0, phiBin(hi - 2.f * pi));
      } else {
        scanBins(phiBin(lo), phiBin(hi));
      }
    }
  };

  // extend one edge inwards: the pixel nodes that make a matching triplet
  // with it, the best curvature match first; returns the new inner edges
  std::vector<std::pair<float, std::uint32_t>> candidates;
  const auto extend = [&](const std::uint32_t e, std::uint32_t& nLinks,
                          std::vector<std::uint32_t>* newEdges) {
    const std::uint32_t a = edges[e].n1;
    const std::uint32_t b = edges[e].n2;
    const auto [lo, hi] = curvatureRange[e];
    const float tauEdge = tau(positions[a], positions[b]);
    // the barrel parameters for an edge of small |eta|, if set
    const bool barrel =
        m_cfg.pixelBarrelMaxAbsEta > 0.f &&
        std::abs(std::asinh(tauEdge)) < m_cfg.pixelBarrelMaxAbsEta;
    const float tolerance = (barrel ? m_cfg.pixelMaxDeltaCurvatureBarrel
                                    : m_cfg.pixelMaxDeltaCurvature) +
                            (barrel ? m_cfg.pixelMaxRelCurvatureBarrel
                                    : m_cfg.pixelMaxRelCurvature) *
                                std::max(std::abs(lo), std::abs(hi));
    const float tauTolerance =
        (barrel ? m_cfg.pixelMaxDeltaTauBarrel : m_cfg.pixelMaxDeltaTau) +
        (barrel ? m_cfg.pixelMaxRelDeltaTauBarrel : m_cfg.pixelMaxRelDeltaTau) *
            std::abs(tauEdge);
    const float phiWindow =
        barrel ? m_cfg.pixelPhiWindowBarrel : m_cfg.pixelPhiWindow;
    const std::uint32_t maxPerEdge =
        barrel ? m_cfg.pixelMaxPerEdgeBarrel : m_cfg.pixelMaxPerEdge;
    candidates.clear();
    forPixelsNear(a, tauEdge - tauTolerance, tauEdge + tauTolerance, phiWindow,
                  [&](const std::uint32_t p) {
                    if (positions[a][3] - positions[p][3] < minDeltaR) {
                      return;
                    }
                    if (std::abs(tau(positions[p], positions[a]) - tauEdge) >
                        tauTolerance) {
                      return;
                    }
                    const float k =
                        curvature(positions[p], positions[a], positions[b]);
                    if (k < lo - tolerance || k > hi + tolerance) {
                      return;
                    }
                    candidates.emplace_back(std::abs(k - 0.5f * (lo + hi)), p);
                  });
    std::ranges::sort(candidates);
    if (candidates.size() > maxPerEdge) {
      candidates.resize(maxPerEdge);
    }
    for (const auto& [score, p] : candidates) {
      const std::uint32_t inner = edgeOf(p, a);
      rawLinks.emplace_back(e, score, inner);
      widen(inner, curvature(positions[p], positions[a], positions[b]));
      ++nLinks;
      if (newEdges != nullptr) {
        newEdges->push_back(inner);
      }
    }
  };

  // strip edges that start a triplet, then the pixel to strip edges made
  const std::uint32_t nStripEdges = static_cast<std::uint32_t>(edges.size());
  std::vector<std::uint32_t> pixelStripEdges;
  for (std::uint32_t e = 0; e < nStripEdges; ++e) {
    if (std::isnan(curvatureRange[e].first) || isPixel(edges[e].n1)) {
      continue;
    }
    if (linksFromTable ? !reachesPixels[nodeView.layers[edges[e].n1]]
                       : positions[edges[e].n1][3] > m_cfg.pixelMaxStripR) {
      continue;
    }
    extend(e, stats.pixelStripLinks, &pixelStripEdges);
  }
  std::ranges::sort(pixelStripEdges);
  const auto [first, last] = std::ranges::unique(pixelStripEdges);
  pixelStripEdges.erase(first, last);
  for (const std::uint32_t e : pixelStripEdges) {
    extend(e, stats.pixelPixelLinks, nullptr);
  }
}

detail::GbtsGraph<detail::TripletGbtsEdge>
TripletGbtsGraphBuilder::buildTheGraph(const TripletGbtsNodes& nodes,
                                       const float bFieldInZ) const {
  constexpr std::uint32_t kNoNode = std::numeric_limits<std::uint32_t>::max();

  Stats stats;

  detail::GbtsGraph<EdgeType> graph;

  const std::uint32_t nNodes = nodes.size();
  const detail::GbtsNodeView nodeView = nodes.nodeView();

  // the strip nodes enter the search in node order: a grid bin sorts its space
  // points by radius, and space points of equal radius then keep this order
  std::vector<std::uint32_t> searchedNodes;
  searchedNodes.reserve(nNodes);
  for (std::uint32_t n = 0; n < nNodes; ++n) {
    if (nodeView.strip(n) != nullptr) {
      searchedNodes.push_back(n);
    }
  }

  // one edge per distinct pair of nodes, inner node first, numbered in the
  // order the pairs first appear in the triplets (inner pair, then outer pair
  // of each triplet). Numbered by sorting rather than with a hash map, which
  // would allocate one node per edge.
  std::vector<EdgeType>& edges = graph.edgeStorage;
  const auto keyOf = [](const std::uint32_t inner, const std::uint32_t outer) {
    return (static_cast<std::uint64_t>(inner) << 32) |
           static_cast<std::uint64_t>(outer);
  };
  const auto makeEdge = [&](const std::uint32_t inner,
                            const std::uint32_t outer) {
    const auto& p1 = nodeView.positions[inner];
    const auto& p2 = nodeView.positions[outer];
    const float dr = p2[3] - p1[3];
    const float tau = dr != 0.0f ? (p2[2] - p1[2]) / dr : 0.0f;
    edges.emplace_back(inner, outer, std::hypot(1.0f, tau) - tau);
  };

  // the layer pairs of the connection table, outer layer first
  std::vector<std::uint8_t> allowedPairs;
  std::uint32_t nLayers = 0;
  if (m_cfg.stripLayersFromTable || m_cfg.stripSearchFromTable ||
      (m_cfg.pixelExtension &&
       (m_cfg.pixelLayersFromTable || m_cfg.pixelLinksFromTable))) {
    const GbtsGeometry& geometry = *m_geometry;
    nLayers = static_cast<std::uint32_t>(geometry.numLayers());
    allowedPairs.assign(static_cast<std::size_t>(nLayers) * nLayers, 0);
    for (const GbtsLayerConnection& c : geometry.layerConnections()) {
      const auto src = geometry.layerIndex(c.src);
      const auto dst = geometry.layerIndex(c.dst);
      if (src.has_value() && dst.has_value()) {
        allowedPairs[static_cast<std::size_t>(*src) * nLayers + *dst] = 1;
      }
    }
  }

  // the triplets, as nodes
  std::vector<std::array<std::uint32_t, 3>> nodeTriplets;
  std::vector<float> nodeTripletQuality;
  m_search->findTriplets(
      {.positions = nodeView.positions,
       .variances = nodes.variances(),
       .strips = nodes.rawStrips(),
       .stripIndex = nodes.stripIndices(),
       .nodes = searchedNodes,
       .offset = m_offset,
       .layers = nodeView.layers,
       .allowedPairs = m_cfg.stripSearchFromTable
                           ? std::span<const std::uint8_t>(allowedPairs)
                           : std::span<const std::uint8_t>(),
       .nLayers = nLayers},
      bFieldInZ, nodeTriplets, nodeTripletQuality, &stats.tripletsBelowQuality);
  stats.triplets = static_cast<std::uint32_t>(nodeTriplets.size());

  // only the triplets whose inner and outer pairs the table connects
  if (m_cfg.stripLayersFromTable) {
    const auto connected = [&](const std::uint32_t inner,
                               const std::uint32_t outer) {
      return allowedPairs[static_cast<std::size_t>(nodeView.layers[outer]) *
                              nLayers +
                          static_cast<std::size_t>(nodeView.layers[inner])] !=
             0;
    };
    std::size_t kept = 0;
    for (std::size_t t = 0; t < nodeTriplets.size(); ++t) {
      const auto& [n0, n1, n2] = nodeTriplets[t];
      if (connected(n0, n1) && connected(n1, n2)) {
        nodeTriplets[kept] = nodeTriplets[t];
        nodeTripletQuality[kept] = nodeTripletQuality[t];
        ++kept;
      }
    }
    stats.tripletsOutsideTable =
        static_cast<std::uint32_t>(nodeTriplets.size() - kept);
    nodeTriplets.resize(kept);
    nodeTripletQuality.resize(kept);
  }

  // occurrence 2t is the inner pair of triplet t, 2t + 1 its outer pair
  const auto nOccurrences = static_cast<std::uint32_t>(2 * nodeTriplets.size());
  std::vector<std::pair<std::uint64_t, std::uint32_t>> occurrences(
      nOccurrences);
  for (std::uint32_t t = 0; t < nodeTriplets.size(); ++t) {
    const auto& [n0, n1, n2] = nodeTriplets[t];
    occurrences[2 * t] = {keyOf(n0, n1), 2 * t};
    occurrences[2 * t + 1] = {keyOf(n1, n2), 2 * t + 1};
  }
  std::ranges::sort(occurrences);
  // the distinct pairs in key order, with the occurrence they first appear at
  std::vector<std::uint64_t> sortedKeys;
  std::vector<std::pair<std::uint32_t, std::uint32_t>>
      firstSeen;  // (occurrence, key rank)
  sortedKeys.reserve(nOccurrences);
  firstSeen.reserve(nOccurrences);
  for (std::uint32_t i = 0; i < nOccurrences; ++i) {
    if (i == 0 || occurrences[i].first != occurrences[i - 1].first) {
      firstSeen.emplace_back(occurrences[i].second,
                             static_cast<std::uint32_t>(sortedKeys.size()));
      sortedKeys.push_back(occurrences[i].first);
    }
  }
  std::ranges::sort(firstSeen);
  std::vector<std::uint32_t> edgeOfKeyRank(sortedKeys.size());
  edges.reserve(sortedKeys.size() + nodeTriplets.size() / 4);
  for (std::uint32_t id = 0; id < firstSeen.size(); ++id) {
    const std::uint32_t rank = firstSeen[id].second;
    edgeOfKeyRank[rank] = id;
    const std::uint64_t key = sortedKeys[rank];
    makeEdge(static_cast<std::uint32_t>(key >> 32),
             static_cast<std::uint32_t>(key & 0xffffffffu));
  }
  std::vector<std::uint32_t> edgeOfOccurrence(nOccurrences);
  {
    std::uint32_t rank = 0;
    for (std::uint32_t i = 0; i < nOccurrences; ++i) {
      if (i > 0 && occurrences[i].first != occurrences[i - 1].first) {
        ++rank;
      }
      edgeOfOccurrence[occurrences[i].second] = edgeOfKeyRank[rank];
    }
  }

  // the edges added later (the pixel extension): looked up among the pairs
  // of the triplets first, then among those added
  std::unordered_map<std::uint64_t, std::uint32_t> addedEdges;
  const auto edgeOf = [&](const std::uint32_t inner,
                          const std::uint32_t outer) {
    const std::uint64_t key = keyOf(inner, outer);
    const auto it = std::ranges::lower_bound(sortedKeys, key);
    if (it != sortedKeys.end() && *it == key) {
      return edgeOfKeyRank[static_cast<std::size_t>(it - sortedKeys.begin())];
    }
    const auto [added, inserted] =
        addedEdges.try_emplace(key, static_cast<std::uint32_t>(edges.size()));
    if (inserted) {
      makeEdge(inner, outer);
    }
    return added->second;
  };

  // the links, outer edge first, with the quality of their triplet
  std::vector<std::tuple<std::uint32_t, float, std::uint32_t>> rawLinks;
  rawLinks.reserve(nodeTriplets.size() + nodeTriplets.size() / 4);
  for (std::uint32_t t = 0; t < nodeTriplets.size(); ++t) {
    rawLinks.emplace_back(edgeOfOccurrence[2 * t + 1], -nodeTripletQuality[t],
                          edgeOfOccurrence[2 * t]);
  }

  if (m_cfg.pixelExtension) {
    // the pixel extension follows the table only with pixelLayersFromTable
    // or pixelLinksFromTable
    const bool pixelTable =
        m_cfg.pixelLayersFromTable || m_cfg.pixelLinksFromTable;
    extendIntoPixels(nodeView, nNodes, edges, edgeOf, rawLinks,
                     pixelTable ? std::span<const std::uint8_t>(allowedPairs)
                                : std::span<const std::uint8_t>(),
                     pixelTable ? nLayers : 0, stats);
  }

  // per outer edge the best triplets first; a link made twice counts once,
  // and an edge links to at most `maxInnerNeighbours` inner edges
  std::ranges::sort(rawLinks);
  std::vector<std::pair<std::uint32_t, std::uint32_t>> links;
  links.reserve(rawLinks.size());
  // the quality of the grid triplet each link was made from; a pixel extension
  // link (its inner edge starts at a pixel node) was made from none
  std::vector<float> linkQuality;
  linkQuality.reserve(rawLinks.size());
  std::vector<std::uint32_t> taken;
  std::uint32_t currentOuter = kNoNode;
  for (const auto& [outer, negQuality, inner] : rawLinks) {
    if (outer != currentOuter) {
      currentOuter = outer;
      taken.clear();
    }
    if (std::ranges::find(taken, inner) != taken.end()) {
      continue;
    }
    if (m_cfg.maxInnerNeighbours != 0 &&
        taken.size() >= m_cfg.maxInnerNeighbours) {
      ++stats.linksOverCap;
      continue;
    }
    taken.push_back(inner);
    links.emplace_back(outer, inner);
    linkQuality.push_back(nodeView.strip(edges[inner].n1) == nullptr
                              ? std::numeric_limits<float>::quiet_NaN()
                              : -negQuality);
  }

  graph.nDoublets = static_cast<std::uint32_t>(edges.size());
  stats.edgesLeftOut = detail::orderOutsideIn(graph, links, linkQuality);
  stats.edges = static_cast<std::uint32_t>(graph.edgeStorage.size());
  stats.links = static_cast<std::uint32_t>(graph.neighbours.size());

  if (m_stats != nullptr) {
    *m_stats = stats;
  }
  return graph;
}

namespace detail {

namespace {
/// An edge not given a place yet
constexpr std::uint32_t kNoEdge = std::numeric_limits<std::uint32_t>::max();
}  // namespace

/// Hand the edges and their links to the graph with every edge before the
/// edges it links to, the order the seeder reads it in.
///
/// An edge is created the first time a triplet passes through it, so an outer
/// edge can come after an inner edge an earlier triplet created. The edges are
/// renumbered, keeping the order they were created in as far as the links
/// allow: unchanged where it already runs outside in. Nor does the cut on the
/// chord rule out a closed loop of links: the edges of a loop, and those only
/// reached through one, are left out with their links.
/// @param graph The graph, its edges in the order they were created
/// @param links The links as they were made, outer edge first
/// @param linkQuality Quality of each link, in the order of `links`; empty for
///        none
/// @return The number of edges left out
std::uint32_t orderOutsideIn(
    detail::GbtsGraph<detail::TripletGbtsEdge>& graph,
    const std::vector<std::pair<std::uint32_t, std::uint32_t>>& links,
    std::span<const float> linkQuality) {
  std::vector<detail::TripletGbtsEdge>& edges = graph.edgeStorage;
  const auto nEdges = static_cast<std::uint32_t>(edges.size());

  // the inner neighbours of each edge, in the order they were linked
  std::vector<std::uint32_t> offsets(nEdges + 1, 0);
  std::vector<std::uint32_t> nOuter(nEdges, 0);
  for (const auto& [outer, inner] : links) {
    ++offsets[outer + 1];
    ++nOuter[inner];
  }
  for (std::uint32_t e = 0; e < nEdges; ++e) {
    offsets[e + 1] += offsets[e];
  }
  std::vector<std::uint32_t> inners(links.size());
  const bool withQuality = !linkQuality.empty();
  std::vector<float> qualities(withQuality ? links.size() : 0);
  std::vector<std::uint32_t> nextSlot(offsets.begin(), offsets.end() - 1);
  for (std::size_t l = 0; l < links.size(); ++l) {
    const auto& [outer, inner] = links[l];
    if (withQuality) {
      qualities[nextSlot[outer]] = linkQuality[l];
    }
    inners[nextSlot[outer]++] = inner;
  }

  // Already outside in, as it mostly is: the order below would leave every
  // edge where it is, so the links are handed over as they are.
  if (std::ranges::all_of(
          links, [](const auto& link) { return link.first < link.second; })) {
    graph.neighbourOffsets = std::move(offsets);
    graph.neighbours = std::move(inners);
    graph.linkQuality = std::move(qualities);
    return 0;
  }

  // Place the earliest created edge whose outer edges are all placed.
  std::priority_queue<std::uint32_t, std::vector<std::uint32_t>, std::greater<>>
      ready;
  for (std::uint32_t e = 0; e < nEdges; ++e) {
    if (nOuter[e] == 0) {
      ready.push(e);
    }
  }
  std::vector<std::uint32_t> order;
  order.reserve(nEdges);
  std::vector<std::uint32_t> newIndex(nEdges, kNoEdge);
  while (!ready.empty()) {
    const std::uint32_t e = ready.top();
    ready.pop();
    newIndex[e] = static_cast<std::uint32_t>(order.size());
    order.push_back(e);
    for (std::uint32_t k = offsets[e]; k < offsets[e + 1]; ++k) {
      if (--nOuter[inners[k]] == 0) {
        ready.push(inners[k]);
      }
    }
  }

  std::vector<detail::TripletGbtsEdge> ordered;
  ordered.reserve(order.size());
  graph.neighbourOffsets.assign(order.size() + 1, 0);
  graph.neighbours.clear();
  graph.neighbours.reserve(links.size());
  graph.linkQuality.clear();
  graph.linkQuality.reserve(qualities.size());
  for (std::uint32_t i = 0; i < order.size(); ++i) {
    const std::uint32_t e = order[i];
    ordered.push_back(edges[e]);
    for (std::uint32_t k = offsets[e]; k < offsets[e + 1]; ++k) {
      // an inner neighbour in a loop was never placed
      if (newIndex[inners[k]] != kNoEdge) {
        graph.neighbours.push_back(newIndex[inners[k]]);
        if (withQuality) {
          graph.linkQuality.push_back(qualities[k]);
        }
      }
    }
    graph.neighbourOffsets[i + 1] =
        static_cast<std::uint32_t>(graph.neighbours.size());
  }
  edges = std::move(ordered);

  return nEdges - static_cast<std::uint32_t>(order.size());
}

}  // namespace detail

}  // namespace Acts::Experimental
