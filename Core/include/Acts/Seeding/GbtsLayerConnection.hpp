// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/GbtsLayerDescription.hpp"

#include <span>

namespace Acts::Experimental {

/// A pair of layers the seeder may build a graph edge between. Edges are made
/// outside-in, so a hit on @c src is the outer end of the doublet and a hit on
/// @c dst the inner one.
///
/// A triplet is cut at the larger of the tolerances its inner and outer pair
/// carry for that role.
struct GbtsLayerConnection {
  /// Outer layer id.
  GbtsExperimentLayerId src{};
  /// Inner layer id.
  GbtsExperimentLayerId dst{};
  /// Largest |tau ratio| if this is the inner pair of a triplet.
  float tauRatioCutInner = 0.007f;
  /// Largest |tau ratio| if this is the outer pair of a triplet.
  float tauRatioCutOuter = 0.007f;
  /// Added to the displaced graph's tau ratio cut for a triplet with a strip
  /// node, by the larger of what its two pairs carry, until the triplet has
  /// resolved its strips.
  float tauRatioCorrStrip = 0.f;
  /// The same, in place of `tauRatioCorrStrip`, between two doublets whose taus
  /// both have their strip ends settled.
  float tauRatioCorrStripSettled = 0.f;
};

/// Tau ratio tolerances for `assignAdaptiveTauRatioCuts`.
struct GbtsAdaptiveTauRatioCuts {
  /// Tolerance of a triplet nothing widens.
  float tauRatioCut = 0.007f;
  /// Widening for a barrel layer skipped, or for the barrel left.
  float tauRatioCorr = 0.006f;
  /// Widening for a strip layer in the triplet.
  float tauRatioCorrStrip = 0.03f;
  /// Widening for a strip layer in the triplet between two doublets whose taus
  /// both have their strip ends settled.
  ///
  /// @note Displaced only.
  float tauRatioCorrStripSettled = 0.0075f;
  /// Fill the pairs for the displaced graph. Every barrel layer then counts as
  /// barrel for `tauRatioCorr`, rather than only the pixel ones, so a pair of
  /// strip barrel layers that skips one is widened, and going from the pixel
  /// barrel onto a strip barrel layer is not leaving it. The strip widenings
  /// go in the pair's own fields rather than in its cuts, to come off once a
  /// triplet has resolved its strips.
  bool displaced = false;
};

/// Fill the tau ratio tolerances of the layer pairs. A pair is widened by
/// `tauRatioCorr` if it skips a barrel layer of any technology between two
/// barrel layers, or, as an outer pair, if it leaves the barrel, the barrel
/// being the pixel one unless `displaced` is set. A pair with a strip layer is
/// widened by `tauRatioCorrStrip` on top, or with `displaced` carries it and
/// `tauRatioCorrStripSettled` in its own fields.
///
/// @param connections The layer pairs, whose tolerances are overwritten
/// @param layers The layers the geometry is built from
/// @param cuts The tolerances to assign
void assignAdaptiveTauRatioCuts(std::span<GbtsLayerConnection> connections,
                                std::span<const GbtsLayerDescription> layers,
                                const GbtsAdaptiveTauRatioCuts& cuts = {});

}  // namespace Acts::Experimental
