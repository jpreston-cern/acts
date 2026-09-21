@defgroup gbts Graph-Based Track Seeding
@ingroup seeding
@brief Seeding by building and filtering a graph of hit doublets

> [!tip]
> This page documents @ref Acts::Experimental::GraphBasedTrackSeeder "GBTS" as
> implemented in ACTS today. GBTS is an alternative to the classical
> triplet-based @ref seeding, not a layer on top of it — the two are independent
> entry points producing the same @ref Acts::SeedContainer.

## Why a graph?

Classical ACTS seeding (@ref seeding) enumerates *triplets* of space points from
a binned grid and cuts on the helix they describe. The combinatorics of that
enumeration grow steeply with occupancy, and every triplet is judged in
isolation.

GBTS inverts the order. It first builds a *graph*: nodes are space points, and a
directed edge joins two space points on connected detector layers whenever the
pair passes a set of cheap two-point cuts. Only then does it look for structure
in that graph — long chains of mutually compatible edges — and turn the best
chains into seeds. The expensive per-candidate work is therefore done once per
*edge* rather than once per triplet, and the chain length itself becomes a
quality signal.

The workflow has four stages, each documented below:

1. @ref gbts-nodes — sort the space points into eta/phi-ordered graph nodes.
2. @ref gbts-graph — create edges between compatible node pairs and link
   compatible edges to each other. @ref gbts-lrt does the same without assuming
   the track came from the beamline.
3. @ref gbts-cca — propagate a "level" through the edge graph to find the
   longest chains.
4. @ref gbts-extraction — follow the best chains through a Kalman-like filter
   and emit seeds.

## Geometry and layer connections {#gbts-geometry}

GBTS does not use the ACTS tracking geometry. It works on its own lightweight
description: a flat list of `GbtsLayer` logical layers, each subdivided into
**eta bins**.

@ref Acts::Experimental::GbtsLayerDescription gives a layer its ID, its type
(barrel or endcap), its sensor technology and its extent. For a barrel layer
`refCoord` is the radius and the bounds are in @f$z@f$; for an endcap it is the
other way round. The ID is the caller's own numbering, and the algorithm never
decodes it.

A barrel layer carries one more field, `depth`: how deep it sits, counting
outwards from zero at the innermost barrel layer the geometry was given. An
endcap layer keeps the default `-1`, so the sign of the field also says whether
a layer is in the barrel. If the caller leaves it unset,
@ref Acts::Experimental::GbtsGeometry sorts the barrel layers by `refCoord` and
fills it in. Set it on every barrel layer or on none of them, because the
constructor rejects a partial ordering.

Three things need to know where a layer sits: the adaptive @f$\tau@f$
correction of @ref gbts-graph, which asks whether three layers are radially
consecutive, and the two innermost-layer cuts of the same section, which ask how
deep a layer is. All three read `depth`, and where a cut means *pixels*
specifically it asks `technology` as well. GBTS therefore runs on any layer
numbering, and on a detector of any technology: `depth` counts every barrel
layer, so it means the same thing whether the innermost layers are pixels or
strips.

> [!note]
> One reader of the ATLAS numbering survives, and it sits outside the core
> algorithm: the examples algorithm decodes the volume id, to pick the strip
> layers out of an ATLAS connection table.

Which layer pairs may be joined by an edge is a list of
@ref Acts::Experimental::GbtsLayerConnection, each naming a source (outer) and a
destination (inner) layer. @ref Acts::Experimental::GbtsGeometry combines the
layer descriptions with those connections and precomputes, for every pair of
connected layers, which *eta bin* pairs are geometrically compatible with the
allowed @f$z_0@f$ range. The result is a **bin group** list — one inner bin
together with all outer bins it may connect to — which serves as the graph
builder's iteration schedule, ordered so that outer bins are processed before
the inner bins that depend on them.

The binning it worked out is readable back off the geometry, so a consumer that
runs the same algorithm elsewhere does not have to recompute or pre-generate it:
@ref Acts::Experimental::GbtsGeometry::layerBinning gives a layer's
@ref Acts::Experimental::GbtsLayerBinning, the eta bins it owns in the global
numbering, and @ref Acts::Experimental::GbtsGeometry::binGroups the schedule
itself. Eta bins are numbered globally and a layer's are contiguous, so the
layers tile the numbering in order. This is what the GPU implementation is
configured from.

> [!note]
> The connections are trained offline rather than written by hand.
> @ref Acts::Experimental::GbtsLayerConnectionTool accumulates layer-pair
> statistics from simulated tracks; the
> `Examples/Scripts/Python/gbts_layer_connection_training_itk.py` and
> `gbts_layer_connection_training_odd.py` scripts drive it for the ITk and the
> Open Data Detector. `ActsExamples::GraphBasedSeedingAlgorithm` reads the
> resulting table, in ATLAS' connector file format, and hands the pairs it
> lists to the geometry.

## Graph nodes {#gbts-nodes}

@ref Acts::Experimental::GbtsNodeStorage holds the graph nodes. Space points are
fed in one at a time through `insert`, which takes **plain scalars** rather than
an ACTS container:

@snippet{trimleft} include/Acts/Seeding/GbtsNodeStorage.hpp gbts insert

As with @ref Acts::CylindricalSpacePointGrid, an experiment can therefore fill
the storage straight from its own space point EDM. Overloads exist for callers
that already have @f$r@f$ and @f$\phi@f$, and for an
@ref Acts::ConstSpacePointProxy together with the columns carrying the layer
index, cluster width and local @f$y@f$ position.

Two different layer numbers meet here, and they are separate types:

| Type | What it is |
| --- | --- |
| `GbtsExperimentLayerId` | the layer id the experiment assign. Sparse and structured -- the layer descriptions and connections are written in terms of it. The algorithm treats it as an opaque key. |
| `GbtsLayerIndex` | where that layer sits in one @ref Acts::Experimental::GbtsGeometry, dense from zero. It indexes the geometry, and it is what a node carries. |

`insert` takes the **index**, because it is on the per-space-point path.
@ref Acts::Experimental::GbtsGeometry::layerIndex hands it out. It is the
geometry's own numbering and is not derivable from the id.

`insert` assigns the node to an eta bin via `GbtsLayer::getEtaBin` and buffers
it. `finalize` then sorts each bin by @f$\phi@f$ and materialises the nodes into
a space point container ordered by eta bin and then by @f$\phi@f$, so that every
eta bin is **one contiguous range of node indices**. A node index is therefore
all that the rest of the algorithm needs to pass around.

The per-node data the graph builder reads lives in dynamic columns on that same
container. It is packed rather than split into one array per field, because the
innermost loop reads all of it together:

@snippet{trimleft} include/Acts/Seeding/detail/GbtsGraphTypes.hpp gbts node params

@f$\tau = \cot\theta@f$. The infinite defaults disable the @f$\tau@f$ cut
entirely; only the optional machine-learning lookup table narrows them (see
@ref gbts-ml). Alongside it sits the bookkeeping the graph builder writes:

@snippet{trimleft} include/Acts/Seeding/detail/GbtsGraphTypes.hpp gbts node edge info

Each eta bin carries its node range plus the @f$\phi@f$ index used by the sliding
window. The @f$\phi@f$ index duplicates entries shifted by @f$\pm 2\pi@f$ near the
wrap-around, so the window never has to handle wrapping:

@snippet{trimleft} include/Acts/Seeding/detail/GbtsGraphTypes.hpp gbts eta bin info

## Building the graph {#gbts-graph}

The builder walks the bin groups from @ref gbts-geometry. For each inner bin it
prepares one **sliding window** in @f$\phi@f$ per connected outer bin, whose
half-width grows with the radial separation of the two bins — the further apart
they are, the more a low-@f$p_T@f$ track can bend between them. It then loops
over the inner nodes, and for each one scans only the outer nodes inside the
window.

A candidate pair @f$(n_1, n_2)@f$ becomes an edge if it survives, in order:

| Cut | Meaning |
| --- | --- |
| @f$\Delta r > @f$ `minDeltaRadius` | the two hits are radially separated enough for @f$\tau@f$ to be meaningful |
| @f$\lvert\tau\rvert < @f$ `maxAbsTau` | within the detector's angular acceptance |
| @f$\tau@f$ inside both nodes' windows | per-node acceptance from the ML lookup table |
| @f$z_0@f$ inside `[minZ0, maxZ0]`, and @f$z@f$ at the outer radius inside the ROI | the pair points back to the luminous region |
| @f$\lvert\kappa\rvert@f$ below an @f$\eta@f$-dependent bound | consistent with the @f$p_T@f$ threshold |

with @f$\tau = \Delta z/\Delta r@f$, @f$z_0 = z_1 - r_1\tau@f$ and the curvature
proxy @f$\kappa = (\phi_2-\phi_1)/\Delta r@f$.

Surviving pairs are appended to a flat edge array:

@snippet{trimleft} include/Acts/Seeding/detail/GbtsGraphTypes.hpp gbts edge

The three fit parameters `p` are @f$\{\exp(-\eta),\ \kappa,\ \phi_1 + \kappa
r_1\}@f$.

Because the inner node's edges are written contiguously, the edges *incoming* to
a node form a contiguous range, recorded in that node's `GbtsNodeEdgeInfo`.
Immediately after creating an edge @f$(n_1, n_2)@f$, the builder scans the edges
incoming to @f$n_2@f$ — that is, edges @f$(n_2, n_3)@f$ — and links the two
whenever the implied triplet is consistent: the @f$\tau@f$ ratio, the @f$\phi@f$
continuation and the curvature difference must all agree within tolerance. For
pixel-barrel triplets an optional
@ref Acts::Experimental::GraphBasedTrackSeeder "validateTriplets" step also fits
a circle through the three points and cuts on @f$d_0@f$ and @f$p_T@f$. Each
edge stores up to `kGbtsMaxEdgeNeighbours` (6) such neighbours.

When `useAdaptiveCuts` is on, the @f$\tau@f$ tolerance is widened for a triplet
that skipped a layer, there being more material to scatter in. Adjacency is
counted over `depth`, so a strip layer sitting between two pixel layers counts
as the gap it is; a triplet whose outermost node left the barrel cannot be
counted at all and is widened on that basis.

Two further cuts apply on the innermost barrel layers, where the combinatorics
are worst. Each has its own limit on the `depth` of the inner layer, and a
negative limit switches the cut off:

- `matchBeforeCreate` (off by default, limited by `matchBeforeCreateMaxDepth`)
  demands the @f$\tau@f$ half of the triplet test *before* the edge exists:
  @f$n_2@f$ must already carry an incoming edge whose @f$\tau@f$ agrees with
  the candidate's within `tauRatioPrecut`. A node with `matchBeforeCreateMaxEdges`
  or fewer incoming edges passes unconditionally, there being too little
  evidence to reject it, so the cut only bites where a node is genuinely busy.
- Every inner node accumulates a 16-bit @f$z_0@f$ **histogram bitmask** of its
  confirmed edges. On the layers down to `z0HistogramMaxDepth` that mask
  rejects candidates whose @f$z_0@f$ falls in an empty bin, and nodes with no
  connections at all are skipped outright. This one is prompt-only: a displaced
  track has no @f$z_0@f$ worth histogramming.

## Connected component analysis {#gbts-cca}

With the edge graph built, a cellular automaton assigns each edge a **level**:
the length of the longest chain of linked edges ending at it. All edges start at
level 1; in each iteration an edge whose level equals that of one of its
neighbours proposes an increment, and the proposals are committed at the end of
the iteration. The sweep repeats until nothing changes, or for at most 15
iterations.

The level is the chain-length signal that drives extraction: an edge at level
@f$L@f$ is the head of a chain spanning @f$L+1@f$ space points.

## Seed extraction {#gbts-extraction}

Edges whose level clears the minimum chain length become **chain heads**, sorted
by level so the longest chains are collected first. Each head is then followed
back through the graph by @ref Acts::Experimental::GbtsTrackingFilter.

The filter is a small Kalman filter over the chain. It carries a state of two
independent parts — a quadratic in the bending plane and a linear @f$z@f$ versus
@f$r@f$ model — and at each step extrapolates to the next node, forms a
@f$\chi^2@f$ residual for each part and rejects the branch if either exceeds its
threshold (`maxDChi2X`, `maxDChi2Y`).

Every accepted hit adds a fixed reward `addHit` to the branch score, minus its
two @f$\chi^2@f$ increments weighted by `weightX` and `weightY`. The score
therefore counts the hits on the chain, discounted by how badly they fit the
circle and the @f$z@f$ versus @f$r@f$ line. Where an edge has several
neighbours the filter *branches*, recursing into each; the branch with the best
accumulated score wins.

The result is a set of seed candidates. These are reduced in two passes:

- **Clone removal.** Candidates are ranked by quality, and each space point is
  assigned to the best candidate claiming it. A candidate that has lost more than
  `hitShareThreshold` of its hits to better candidates is dropped.
- **Seed splitting.** Short, central candidates are checked for self-consistency
  by fitting the circle through three different hit subsets. If the three
  curvature estimates disagree by more than `maxInvRadDiff`, the candidate is
  emitted as two shorter "drop-out" seeds instead of one.

The surviving candidates are written to the output @ref Acts::SeedContainer, with
node indices translated back to the caller's own space point indices.

## Machine-learning assisted acceptance {#gbts-ml}

When `useClusterWidthCuts` is enabled, GBTS narrows the per-node @f$\tau@f$
window using a pre-trained lookup table indexed by **pixel cluster width**. The
cluster a track leaves in a pixel module grows with the incidence angle, so the
width alone constrains @f$\cot\theta@f$ before any pairing is attempted.

The table carries two sets of bounds per width bin: one for clusters comfortably
inside the module, and one for clusters within `moduleEdgeTolerance` of the module
edge, where the cluster may be truncated and the width therefore underestimates
the angle. Wide clusters in the pixel endcap are dropped entirely
(`maxEndcapClusterWidth`).

> [!note]
> The seeder takes the table itself as `tauLookupTable`, not a path to it;
> `ActsExamples::GraphBasedSeedingAlgorithm` parses it from ATLAS' text format.
> It is only consulted for pixel barrel layers, and the ACTS examples framework
> does not currently provide cluster widths or local positions, so this path is
> exercised only by experiment-side integrations that supply them through
> `insert`.

## Large radius tracking {#gbts-lrt}

Everything above assumes the track came from the luminous region. That
assumption is not a cut that can be loosened; it is load-bearing. The prompt
builder reads a curvature off a *pair* of hits as @f$\kappa =
(\phi_2-\phi_1)/\Delta r@f$, which is only a curvature because the beamline is
silently used as a third point on the circle. Two transverse points do not fix a
circle. Neither do they fix @f$z_0@f$, which the prompt form takes as @f$z_1 -
r_1\tau@f$ — true only for a track leaving the beamline radially.

For a track from a displaced vertex all of that is unavailable at doublet stage.
@ref Acts::Experimental::DisplacedGbtsGraph builds the same kind of graph without
it. It is a drop-in alternative to
@ref Acts::Experimental::GbtsGraphBuilder "GbtsGraphBuilder" — same node
storage, same configuration object, same edge levels and chain heads, so
@ref gbts-cca and @ref gbts-extraction are unchanged and the seeder is written
once over either.

### What a doublet can still say

| Quantity | Prompt | Displaced |
| --- | --- | --- |
| @f$\tau = \cot\theta@f$ | @f$\Delta z/\Delta r@f$ | @f$\Delta z/\mathrm{chord}@f$ |
| curvature | from the pair | **not determined** |
| azimuth at the perigee | from the pair | **not determined** |
| @f$z_0@f$ | a value | **a band** |

So the doublet stage keeps only the @f$r@f$–@f$z@f$ (non-bending) plane, and the
cuts that lived there change shape:

- **The @f$z_0@f$ band.** @f$z@f$ is linear in the transverse path @f$S_1@f$ from
  the perigee to the inner hit, and @f$S_1 = \sqrt{r_1^2 - d_0^2}@f$ with
  @f$d_0@f$ unknown. @f$S_1@f$ is monotonic in @f$d_0@f$, so the two ends of
  @f$d_0 \in [0, d_0^{\max}]@f$ bracket a band. The band has to *overlap*
  @f$[z_{0}^{\min}, z_{0}^{\max}]@f$ rather than sit inside it. The same holds
  for @f$z@f$ at the outer radius.
- **The @f$\phi@f$ window.** A track of impact parameter @f$d_0@f$ sits at
  azimuth @f$\arccos(d_0/r)@f$, so between two radii it walks
  @f$\arccos(d_0^{\max}/r_2) - \arccos(d_0^{\max}/r_1)@f$ on displacement alone,
  on top of whatever it bends. Where the track started and which way it bends are
  independent, so the window is the **sum** of the two terms, not their
  difference. A layer inside the @f$d_0@f$ limit takes the limit down to its own
  radius, since a track only reaches radius @f$r@f$ if its @f$\lvert d_0\rvert@f$
  is below @f$r@f$.
- **No curvature cut.** There is no curvature yet to cut on.

The price is combinatorics, and `d0Max` sets its scale: it widens the
@f$\phi@f$ window and the @f$z_0@f$ band together.

### The triplet fit {#gbts-lrt-fit}

Three nodes do fix a circle, so the edge-linking step is where the track
parameters first exist. Inverting the transverse plane about the middle node
maps every circle through that node onto a straight line, so the fit is the line
through the two remaining nodes and needs no iteration. It is the same conformal
mapping the prompt graph's `validateTriplets` uses, kept whole rather than
reduced to a verdict: signed curvature, @f$d_0@f$, @f$\cot\theta@f$ along the
fitted arc, and the track tangent at each of the three nodes.

Curvature is reported in the prompt graph's @f$\mathrm{d}\phi/\mathrm{d}r@f$
convention — half the geometric @f$1/R@f$ — so the tuned cut values carry over
unchanged.

**Strip resolution moves here too.** A strip space point has to be slid along its
strip to where the track crossed it, and that needs a direction. The only
direction a doublet has is its chord, which is the tangent *only* for a track
from the beamline. The displaced builder therefore resolves each strip end
against the fitted tangent and refits on the moved ends. One pass is enough: the
slide is a small fraction of a strip and the fit is linear in the node positions
to that order.

The triplet is then cut on what the fit gives: the @f$\tau@f$ ratio retaken on
the resolved points, an @f$\eta@f$-dependent bound on the curvature, and
@f$d_0@f$ and @f$p_T@f$ under `validateTriplets`. The @f$\tau@f$ ratio is taken
twice on purpose — loosely before the fit, keeping `tauRatioCorrStrip`, because
the along-strip coordinate it leans on is the very thing about to move, and
tightly afterwards with that looseness dropped for every end the calibration
actually put back.

### Skipping a candidate before fitting it {#gbts-lrt-turn}

The fit is far too expensive to run on every pair of edges meeting at a node. A
chord bisects the tangents at its ends, so the turn from the inner chord to the
outer one is @f$\arcsin(\kappa L_{13})@f$: it grows with curvature and so falls
with @f$p_T@f$, and a wide turn is not worth fitting.

@f$L_{13}@f$ is not measured — the two chords bound it, @f$L_{13} \le L_{12} +
L_{23}@f$ — and each edge caches its own chord direction as a unit vector when it
is built, so the turn's sine and cosine come from a dot and a cross of two cached
pairs. No trigonometry, no square root, and the angle is never formed. On
adjacent barrel layers it rejects between a third and two thirds of candidates
before the fit sees them.

The cut is skipped where a strip end has still to slide, since in the endcap that
slide is largely radial and moves the very chords it reads. Rejecting a candidate
the fit would have kept is the one thing it must not do.

### Matching triplets {#gbts-lrt-matching}

The prompt graph links two doublets by comparing the parameters each one carries.
A displaced doublet carries none, so an edge instead records the triplets it is
the **inner** edge of:

@snippet{trimleft} include/Acts/Seeding/detail/DisplacedGraphTypes.hpp displaced triplet properties

Held inline and counted, exactly as `vNei` is:

@snippet{trimleft} include/Acts/Seeding/detail/DisplacedGraphTypes.hpp displaced edge

Bin groups are walked outer to inner, so when an edge @f$(n_1, n_2)@f$ is made,
every edge hanging off @f$n_2@f$ already exists *and already carries its
properties*. A candidate triplet @f$(n_1, n_2, n_3)@f$ is matched against the
triplets its outer edge @f$(n_2, n_3)@f$ already belongs to. Those share a whole
*doublet* with it, so they are the same track candidate seen one node apart, and
both hold a tangent azimuth at @f$n_2@f$. Agreement in that azimuth, in curvature
and in the arc-corrected @f$\tau@f$ is the displaced stand-in for the prompt
@f$\Delta\phi@f$/@f$\Delta\kappa@f$ test.

The bookkeeping that makes this line up is that an edge's stored @f$\phi@f$ is
*always the tangent at that edge's own inner node*. A triplet writes onto its
inner edge the tangent at node 1, and compares against its outer edge the tangent
at node 2 — the same node, reached from two positions in the triplet, because the
middle node of a triplet is the inner node of its outer edge.

Agreement with **any one** recorded triplet is enough. An outer edge may belong
to several from combinatorics, only one of which is the real track, and the graph
keeps that ambiguity for @ref gbts-cca and @ref gbts-extraction to resolve.

An outer edge with no triplets yet is the outermost pair of its chain and has
nothing to disagree with, so it passes. That is worth knowing when tuning: the
**first triplet of every chain** is accepted on the absolute cuts alone —
@f$d_0@f$, @f$p_T@f$, curvature and the @f$\tau@f$ ratio — with no relative test
behind it.

### Chain length

A link here has already cost a full triplet fit, so two edges are three nodes
that were fitted and cut on rather than three that merely agreed on @f$\tau@f$.
`minSeedLevel` can be set to 2 accordingly, where the prompt graph wants 3.

> [!warning]
> `addTriplets` accepts a chain one level short, so at a `minSeedLevel` of 2 it
> reaches down to a single edge — two nodes and no fit at all. Raise
> `minSeedLevel` before enabling it.

### Running it

`ActsExamples::GraphBasedSeedingAlgorithm` selects the builder with
`useDisplacedGraph`; everything else is shared, including `graphConfig`. Two
settings outside the graph default to a beamspot the displaced mode is meant to
leave behind, and both silently cost tracks if left alone:

- @ref Acts::Experimental::GbtsZ0Range, which fixes which eta bins may be linked
  at all, defaults to @f$\pm 168\,\mathrm{mm}@f$. A track starting further down
  the beamline has no link to travel along however well it passes the cuts.
- @ref Acts::Experimental::GbtsTrackingFilter "GbtsTrackingFilter::Config::maxZ0"
  defaults to @f$170\,\mathrm{mm}@f$ and discards what the graph did find.

## Configuration {#gbts-configuration}

The cuts that build and link the doublets live on
@ref Acts::Experimental::GbtsGraphConfig, which `GraphBasedSeedingAlgorithm`
takes as `graphConfig` and exposes to Python as `GbtsGraphBuilderConfig`.
**Both builders take the same object**, so a caller sets it up once and hands
it to whichever it runs; what differs is not which knobs exist but where each
one bites, and the fields whose meaning moves say so.

Fields marked *prompt* are read by @ref Acts::Experimental::GbtsGraphBuilder
alone.

| Option | Stage | Effect |
| --- | --- | --- |
| `minPt` | @ref gbts-graph | drives the curvature and @f$\phi@f$-window bounds |
| `minDeltaRadius`, `maxAbsTau` | @ref gbts-graph | doublet acceptance |
| `minZ0`, `maxZ0`, `doubletFilterRZ` | @ref gbts-graph | luminous-region cuts on the doublet |
| `tauRatioCut`, `cutDPhiMax`, `cutDCurvMax` | @ref gbts-graph | edge-to-edge linking tolerances |
| `useAdaptiveCuts`, `tauRatioCorr` | @ref gbts-graph | widen the @f$\tau@f$ tolerance when a layer is skipped, counted over `depth` |
| `tauRatioCorrStrip` | @ref gbts-graph | widen it again where an end is an unresolved strip. Displaced, this comes back off once the triplet has resolved it |
| `calibrateStrips`, `maxStripLengthFraction` | @ref gbts-graph | resolve a strip along its strip. Prompt: per doublet, against the chord. Displaced: per triplet, against the fitted tangent |
| `validateTriplets` | @ref gbts-graph | prompt: circle fit on pixel-barrel triplets. Displaced: the @f$d_0@f$ and @f$p_T@f$ cuts of a fit that runs regardless |
| `d0Max` | @ref gbts-graph | largest @f$d_0@f$ a triplet may have. Displaced, it also sets the @f$z_0@f$ band and the @f$\phi@f$ window, so it drives the combinatorics |
| `nMaxEdges` | @ref gbts-graph | hard cap on the edge array (2M by default); exceeding it costs efficiency |
| `matchBeforeCreate`, `tauRatioPrecut`, `matchBeforeCreateMaxDepth`, `matchBeforeCreateMaxEdges` | @ref gbts-graph | require a compatible incoming edge before creating one, down to that depth in the barrel |
| `z0HistogramMaxDepth`, `z0Resolution` | @ref gbts-graph | *prompt*: @f$z_0@f$ histogram cut, down to that depth in the barrel |
| `ccaMaxIterations` | @ref gbts-graph | cap on the connected component iterations |
| `minSeedLevel` | @ref gbts-extraction | chain length a candidate must reach, in edges. 3 by default; a displaced run can use 2, since a link there already cost a triplet fit |
| `addTriplets`, `maxAbsEtaAddTriplets` | @ref gbts-extraction | allow shorter chains within an @f$\eta@f$ range |

The graph applies the last three when it picks the chain heads, so they sit
with the graph rather than with the seeder that reads those chains back out.

The rest are on @ref Acts::Experimental::GraphBasedTrackSeeder "GraphBasedTrackSeeder::Config":

| Option | Stage | Effect |
| --- | --- | --- |
| `nMaxPhiSlice` | @ref gbts-nodes | sets the @f$\phi@f$ slice width, and with it the base sliding-window width the graph uses |
| `hitShareThreshold` | @ref gbts-extraction | fraction of shared hits above which a candidate is a clone |
| `maxSeedSplitEta`, `maxInvRadDiff` | @ref gbts-extraction | seed splitting |
| `useClusterWidthCuts`, `tauLookupTable` | @ref gbts-ml | cluster-width based @f$\tau@f$ windows |
| `maxEndcapClusterWidth`, `moduleHalfLengthY`, `moduleEdgeTolerance` | @ref gbts-ml | cluster-width acceptance and module-edge handling |

`useStripConnections`, which takes the strip layer connections from the
connector file instead of the pixel ones, is read where the file is loaded and
so sits on `GraphBasedSeedingAlgorithm::Config` itself.

@ref Acts::Experimental::GbtsTrackingFilter "GbtsTrackingFilter::Config"
separately controls the chain-following filter of @ref gbts-extraction "seed extraction":

| Option | Effect |
| --- | --- |
| `sigmaX`, `sigmaY` | measurement resolution in the bending plane and along @f$z@f$ |
| `maxDChi2X`, `maxDChi2Y` | per-step @f$\chi^2@f$ ceilings; a branch exceeding either is dropped |
| `addHit`, `weightX`, `weightY` | the reward and the two @f$\chi^2@f$ weights in the branch score |
| `sigmaMS`, `radLen` | multiple-scattering inflation added before each extrapolation |
| `maxCurvature`, `maxZ0` | track-level bounds checked after each update |

## Implementation pointers {#gbts-implementation}

- Seeder and configuration: @ref Acts::Experimental::GraphBasedTrackSeeder and
  @ref Acts::Experimental::GbtsGraphConfig, shared by both builders. The seeder
  and @ref Acts::Experimental::GbtsTrackingFilter are templated on the graph,
  which names its own edge type, so the call sites deduce it.
- Graph builders: @ref Acts::Experimental::GbtsGraphBuilder and
  @ref Acts::Experimental::DisplacedGbtsGraph. The displaced EDM -
  `DisplacedGbtsEdge`, `TripletProperties` and the `TripletCircle` the fit
  returns - lives in `Acts/Seeding/detail/DisplacedGraphTypes.hpp`.
- Node storage: @ref Acts::Experimental::GbtsNodeStorage. The graph EDM it
  holds - `GbtsNodeParams`, `GbtsNodeEdgeInfo`, `GbtsEtaBinInfo`, `GbtsEdge` -
  is internal and lives in `Acts/Seeding/detail/GbtsGraphTypes.hpp`.
- Geometry: @ref Acts::Experimental::GbtsGeometry,
  @ref Acts::Experimental::GbtsLayerConnection, the binning it hands back in
  `Acts/Seeding/GbtsBinning.hpp`, and the internal `GbtsLayer`.
- Chain following: @ref Acts::Experimental::GbtsTrackingFilter and its internal
  `GbtsEdgeState`.
- Region of interest: @ref Acts::Experimental::GbtsRoiDescriptor.
- Connection-table training: @ref Acts::Experimental::GbtsLayerConnectionTool.
- Examples integration: `ActsExamples::GraphBasedSeedingAlgorithm`, driven from
  `Examples/Scripts/Python/full_chain_itk_Gbts.py`.

A GPU implementation of the same algorithm, using an equivalent
struct-of-arrays layout, lives in the traccc plugin under
`Traccc/device/common/include/traccc/gbts_seeding`.
