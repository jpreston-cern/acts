// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"

#include <cstdint>

namespace Acts::Experimental {

/// Configuration for either GBTS graph builder.
///
/// The prompt and displaced builders take the same settings, so they share one
/// configuration and a caller can hand the same object to whichever it runs.
/// What differs is not which knobs exist but where each one bites: the prompt
/// graph fixes a circle from a doublet by borrowing the beamline as a third
/// point, and the displaced graph, which cannot make that assumption, defers
/// the same quantity to a triplet fit. Fields whose meaning moves that way say
/// so, as do the few one builder ignores.
struct GbtsGraphConfig {
  /// Match seeds before creating them.
  bool matchBeforeCreate = false;

  /// Cut a triplet on the pT and d0 of its own fit. The prompt graph does this
  /// for barrel triplets only; the displaced graph does it for every triplet,
  /// and runs the fit either way since the curvature properties an edge
  /// carries come out of it.
  bool validateTriplets = true;

  /// Widens allowed variation in tau ratio if a layer is missed in edge
  /// connecting.
  bool useAdaptiveCuts = true;

  /// Tau ratio cut threshold.
  float tauRatioCut = 0.007f;

  /// Tau ratio precut threshold.
  float tauRatioPrecut = 0.009f;

  /// Correction applied to tau acceptance if a layer is missed during edge
  /// connecting.
  float tauRatioCorr = 0.006f;

  /// The same for a triplet any of whose three nodes a strip module made,
  /// whose two doublets resolved the shared node's along-strip coordinate
  /// separately. Reaches nothing without a strip in the triplet.
  float tauRatioCorrStrip = 0.03f;

  /// Minimum transverse momentum.
  float minPt = 1.0f * Acts::UnitConstants::GeV;

  /// Fraction of `minPt` a triplet may fall to, allowing for three-point
  /// pT resolution.
  float tripletPtFraction = 0.8f;

  // Graph-building options

  /// Use eta binning from geometry structure.
  bool useEtaBinning = true;

  /// Apply RZ cuts on doublets.
  bool doubletFilterRZ = true;

  /// Maximum number of GBTS edges/doublets.
  std::uint32_t nMaxEdges = 2000000;

  /// Minimum delta radius between layers.
  float minDeltaRadius = 2.0f * Acts::UnitConstants::mm;

  /// Largest |cot(theta)| accepted for a doublet. The default corresponds to
  /// |eta| of about 4.3, beyond the acceptance of any current tracker.
  float maxAbsTau = 36.0f;

  /// Maximum d0 impact parameter when validating an edge-connection triplet.
  ///
  /// @note The displaced graph also reads this at doublet stage, where d0 is
  ///       still unknown: it sets the width of the z0 band and of the phi
  ///       window. Opening it there costs combinatorics, not just acceptance.
  float d0Max = 3.0f * Acts::UnitConstants::mm;

  /// Maximum difference in allowed tangent between candidate edge
  /// connections. Prompt, the two tangents come from the two doublets;
  /// displaced, from the two triplet fits, taken at the doublet they share.
  float cutDPhiMax = 0.012f;

  /// Maximum allowed curvature tolerance for candidate edge connections.
  float cutDCurvMax = 0.001f;

  /// Minimum z0 value. In pixel mode the value is picked from the RoI.
  float minZ0 = -600.0f;

  /// Maximum z0 value. In pixel mode the value is picked from the RoI.
  float maxZ0 = 600.0f;

  /// pT at which the default cut coefficients were tuned; they scale by
  /// `tuningPt / minPt`.
  float tuningPt = 0.9f * Acts::UnitConstants::GeV;

  /// Maximum |curvature| above `curvatureSplitAbsTau`, before that scaling.
  /// Applied to a doublet by the prompt graph and to a fitted triplet by the
  /// displaced one, which has no curvature until then.
  float maxCurvatureHighEta = 4.75e-4f / Acts::UnitConstants::mm;

  /// Maximum |curvature| below `curvatureSplitAbsTau`, before that scaling.
  float maxCurvatureLowEta = 3.75e-4f / Acts::UnitConstants::mm;

  /// |cot(theta)| separating the two curvature cuts, corresponding to |eta|
  /// of about 2.1.
  float curvatureSplitAbsTau = 4.0f;

  /// Radial separation splitting the two phi-window slopes below.
  float phiWindowSplitDeltaRadius = 60.0f * Acts::UnitConstants::mm;

  /// Phi window below `phiWindowSplitDeltaRadius`, as an offset plus a slope
  /// times the radial separation.
  float phiWindowNearOffset = 0.002f;

  /// Slope of the near phi window per unit radial separation. Scaled by
  /// `tuningPt / minPt`.
  float phiWindowNearSlope = 4.33e-4f / Acts::UnitConstants::mm;

  /// Phi window above `phiWindowSplitDeltaRadius`, in the same form.
  float phiWindowFarOffset = 0.015f;

  /// Slope of the far phi window per unit radial separation. Scaled by
  /// `tuningPt / minPt`.
  float phiWindowFarSlope = 2.2e-4f / Acts::UnitConstants::mm;

  /// Incoming edge count below which a node is accepted without a tau match.
  std::uint32_t matchBeforeCreateMaxEdges = 2;

  /// Deepest barrel layer, counting outwards from zero at the innermost the
  /// geometry holds, whose nodes are cut against the z0 histogram of their
  /// outer neighbourhood and whose isolated nodes are skipped. A negative
  /// value disables the cut.
  ///
  /// @note Prompt only. A displaced track has no z0 worth histogramming.
  std::int32_t z0HistogramMaxDepth = 0;

  /// Deepest barrel layer, counting outwards from zero at the innermost the
  /// geometry holds, to which `matchBeforeCreate` applies when it is enabled.
  /// A negative value disables it, and an endcap layer is never included.
  ///
  /// Counted over every barrel layer whatever its technology, so it selects
  /// the same layers on a detector whose innermost are strips as on one whose
  /// innermost are pixels.
  std::int32_t matchBeforeCreateMaxDepth = 1;

  /// Half-width of the z0 window against which a node is matched in the
  /// histogram.
  ///
  /// @note Prompt only, with `z0HistogramMaxBarrelOrder`.
  float z0Resolution = 2.5f * Acts::UnitConstants::mm;

  /// Maximum radius of the pixel detector.
  float maxOuterRadius = 550.0f;

  /// Resolve a strip node along its strip against the track direction before
  /// cutting on it. Nothing is written back; the correction belongs to the
  /// candidate that made it. The prompt graph does this per doublet, against
  /// the chord; the displaced graph per triplet, against the fitted tangent,
  /// the chord being the tangent only for a track from the beamline.
  bool calibrateStrips = true;

  /// How far along a strip a crossing may land and still be recovered, as a
  /// multiple of the strip half-length, so 1 is the strip itself. This is the
  /// same quantity as `TripletSeedFinder::Config::toleranceParam`.
  float maxStripLengthFraction = 1.1f;

  /// Maximum number of connected-component iterations.
  std::uint32_t ccaMaxIterations = 15;

  // Chain selection options, shared with the seed extraction that reads the
  // chains back out of the graph.

  /// Chain length a seed candidate must reach, counted in edges: a triplet
  /// plus one confirmation.
  ///
  /// @note The displaced graph can be run at two, where the prompt graph
  ///       cannot: a link there already cost a triplet fit, so two edges are
  ///       three nodes that were fitted and cut on rather than three that
  ///       merely agreed on tau.
  std::uint32_t minSeedLevel = 3;

  /// optionally add 3 sp seeds within a certain eta range
  ///
  /// @note Worth little until `maxAbsEtaAddTriplets` is opened past
  ///       `edgeMaskMinEta`; matters most where there are few layers.
  /// @note This accepts a chain one level short, so at a `minSeedLevel` of two
  ///       it reaches down to a single edge: two nodes and no fit at all.
  bool addTriplets = false;

  /// the maximum allowed eta value in which
  /// three spacepoint seeds are passed through
  float maxAbsEtaAddTriplets = 1.5;
};

}  // namespace Acts::Experimental
