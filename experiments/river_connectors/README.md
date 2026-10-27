# River connector experiment

Standalone geometry diagnostic, outside the terrain build and water-carver.
Reads cached Jamestown polygons and centrelines; writes only to
`output/jamestown-river-lab/connectors/<from-id>-<to-id>/`.

From the repository root:

```sh
uv run python experiments/river_connectors/run.py --from-id 2240 --to-id 4563
```

Options: `--sample-run 75` metres inside each polygon, `--sample-step 5`
metres, `--tangent-window 10` metres either side of each sample.

The runner selects a single source centreline intersecting both polygons.
The geometry module finds the shortest uncovered along-line gap between them.
It measures perpendicular cross-sections at 5 m intervals, starting 2.5 m
inside each end and continuing up to 75 m (or the available polygon reach).
Only the connected water interval containing the sample point is measured.
Each end's arithmetic mean gets equal weight in the final mean width.

That constant width buffers the original gap route. For overlap at the joins,
the route is extended by up to one connector width inside each polygon,
limited to the adjacent intersected reach. It does not trim, replace or edit
the source polygons. The resulting footprint is a proposal, not a guaranteed
bank-aligned union: offset centrelines and different end widths can leave
shoulders. The original centreline's angular bends are retained.

Outputs: a preview showing every accepted transect, a JSON measurement report,
and a WGS84 GeoJSON containing the proposed corridor, gap line and transects.
Report endpoint means are in route station order; sample `end` values index
`polygon_ids`. No elevations, carving, source downloads or production imports
of this experiment are involved.

Limits: one simple continuous route, shortest gap only, no automatic network
routing or bank-end detection. Curved/oblique sections can bias measured
widths. This first experiment deliberately makes those effects visible.

## Endpoint-aligned smooth curve trial

```sh
uv run python experiments/river_connectors/run_smooth.py
uv run python -m unittest discover -s experiments/river_connectors -p 'test_*.py'
```

This trial is fixed to Awarua 2240/4563. `smooth.fit_curve` is the reusable
geometry function; its distance parameters are in metres. It fits a straight
local trend through transect midpoints within 15 m of an anchor 20 m inside
each polygon. Those trends supply endpoint positions and unit tangents for
a clamped cubic spline. Interior guide points every 40 m follow the original
route after modest 10 m Gaussian smoothing. Output sampling is approximately
1 m in source route station, not exact curve arc length.

The original 32.9 m constant width remains unchanged, to isolate alignment
and curve smoothing from bank blending. Output goes in the separate `smooth/`
subfolder, including a comparison image, geographic geometry and measurements.
The image's top panels share bounds; the bottom panels are equal-scale zooms.
Self-crossing curves are rejected. Neither land avoidance nor maximum offset
is enforced: reported offsets are diagnostics, not solver constraints. The
local fitted anchors are estimates, and sharp bends can bias their directions.
The trial does not claim that buffering the new line perfectly matches banks.

## Addition-only seam smoothing

```sh
uv run python experiments/river_connectors/run_joins.py --radius 5 --seam-reach 25
```

Run after the smooth-curve trial. `joins.smooth_joins` unions the source polygons
and synthetic corridor, then proposes morphological closing with round buffers
(outward by 5 m, inward by 5 m). Only new area within 25 m of the original and
connector boundary intersections is accepted. This also fills nearby concave
source corners inside that zone; protruding corners remain. Locality is a hard
mask, not a gradual fade, so inspect mask-edge seams for other inputs.

The final geometry includes every part of the original polygons AND corridor.
The exported `synthetic-patch` is final minus original, so original plus patch
reconstructs the result. `smoothing-additions` isolates the extra fill for review;
`final-union` is a convenience export. Neither the downloaded sources nor terrain
outputs are changed. Results are in `joins-5m/` (or the chosen radius).
Measurements verify preservation, locality, validity and patch reconstruction,
allowing tiny floating-point overlay tolerances. Previews show before/after on
identical scales, with added fill in orange.

Compare 5, 10, 25 and 50 m radii:

```sh
uv run python experiments/river_connectors/run_join_sweep.py
```

Outputs are separate in `join-radius-sweep/`: equal-scale join panels, a full
connector overview, measurements and a GeoJSON per radius. Seam reach grows
with radius (`max(25, 3 * radius)` metres), so this compares broader local
smoothing treatments, not radius alone with a fixed spatial mask. At 50 m the
zone includes the connector's bend and closing widens that bend substantially.
All variants retain the source and connector areas; no default is changed.

## Next reach: 4563 to 1108

The smooth and join runners now accept polygon IDs and derive their close-up
positions from the fitted endpoints rather than hardcoded map coordinates:

```sh
uv run python experiments/river_connectors/run_smooth.py --from-id 4563 --to-id 1108
uv run python experiments/river_connectors/run_joins.py --from-id 4563 --to-id 1108
```

The join runner's preferred default is now **20 m radius**, with automatic
seam reach `max(25, 3 * radius)` (60 m here); `--seam-reach` overrides it.
The core geometry function retains its explicit/default API independently.
Outputs remain per pair and per radius. `overview.png` shows the whole new
connector; `comparison.png` shows before/after joins at identical scales.
The old statement that the smooth runner is fixed to 2240/4563 is superseded
by these arguments (2240/4563 remain the CLI defaults). No production carving
configuration has changed.

## Widths matched to both ends

```sh
uv run python experiments/river_connectors/run_smooth.py --from-id 4563 --to-id 1108 --width-transition 200
uv run python experiments/river_connectors/run_joins.py --from-id 4563 --to-id 1108 --curve-variant variable-width
```

`widths.tapered_corridor` starts at each end's measured mean width and uses a
smoothstep ramp to their average over 200 m of curve arc length. Ramps start at
the fitted interior anchors (20 m inside the polygons), not the precise boundary.
The middle stays at the average. On short connectors, each ramp is shortened
to half the curve length. Cross-sections are sampled approximately every metre;
the corridor is the union of triangular strips between successive sections,
which resolves overlapping offsets at tight bends like a conventional buffer.
Consequently the nominal width profile is not a guarantee of exact local width
where a bend overlaps itself. Endpoint caps are flat.

These outputs are isolated in `variable-width/` and `variable-width-joins-20m/`;
previous constant-width trials remain available. A width-profile plot documents
the two ramps. Join smoothing still uses 20 m radius and 60 m locality. Source
polygons remain untouched, and seam additions may widen the nominal corridor.

Compare the 200 m transitions against a linear taper along the whole saved curve:

```sh
uv run python experiments/river_connectors/compare_width_modes.py
```

Both variants use the identical 4563/1108 fitted centreline, endpoint widths,
20 m addition-only smoothing and 60 m seam zones. `widths.tapered_corridor`
accepts `mode="linear"` to interpolate by curve arc length rather than using
a middle plateau. Outputs in `width-mode-comparison/` include width profiles,
equal-scale maps and join details, geometry and preservation checks. This is
a comparison only; the existing default mode remains `end-transitions`.

## Current selection and McKenzie inspection

Linear width interpolation is now the experiment default in `tapered_corridor`
and `run_smooth.py`. The latter writes to `linear-width/`; `run_joins.py` defaults
to that variant and writes to `linear-width-joins-20m/`. For historical trials,
use `--width-mode end-transitions --width-transition 200` or
`--width-mode constant`, then select the corresponding join variant explicitly.
This supersedes earlier default descriptions. Free/headwater endpoints (such as
0.5 m starts) are not yet inferred; the width primitive accepts such a positive
width, but the route discovery still requires two polygon intersections.

```sh
uv run python experiments/river_connectors/inspect_mckenzie.py
```

This read-only source diagnostic shows polygons 180, 2238, 2239 and 4562 and
nearby cached LINZ lines, separating named McKenzie features from other lines.
There is no OSM route in the displayed area. Named gap lines meet polygon
boundaries without positive-length interior overlap, unlike the Awarua route.
Two gaps have paired lines meeting opposite bank corners; these source features
are classified river_cl, so interpreting them as banks needs further work.
The existing connector cannot sample inside polygons along these short lines
without an additional endpoint/route construction step. No McKenzie connector
or source/headwater taper has been generated by the inspection.

## McKenzie separate-channel trial

```sh
uv run python experiments/river_connectors/trial_mckenzie.py
```

Trials all five named gap lines independently: 380913/677434 between 180/2238,
677435 between 2238/2239, and 677436/677437 between 2239/4562. Extend each end
75 m along its local 10 m tangent, then reuse width sampling, midpoint alignment,
linear tapering and 20 m addition-only smoothing (60 m seam reach).

Two straight extensions miss their destination polygon. A recorded fallback
rotates only the invented interior extensions in 5-degree increments, selecting
the first direction with 35 m of interior sampling runway. Both need a 5-degree
turn in this dataset. Original gap lines remain intact before curve fitting.
No channel receives a split-flow width reduction: each samples the polygon
independently, as requested. These inferred widths are not observed branch widths.

Outputs in `mckenzie-separate-lines/` contain per-line status, widths, source and
adjusted extensions, samples, patches and their combined footprint. Top panels
show channels separately with transects; bottom panels show the union, including
holes. Smoothing may fill existing land holes within the 60 m zone: the first
gap shows this especially clearly. Addition-only preserves source water, not
source islands. This trial intentionally exposes that limitation for review;
no source files or production heightmaps are modified.

### Accepted McKenzie result

The separate-channel McKenzie trial is the selected approach for these five
lines: retain independently measured branch widths, linear interpolation and
20 m addition-only seam smoothing. The existing generated geometry remains the
accepted review artifact; acceptance does not yet wire it into terrain builds.

Repeatability boundary: geometry operations are computed from input shapes,
not hand-drawn. The McKenzie runner explicitly selects five source line IDs
and their polygon pairs, and its preview windows are location-specific.
Straight-extension fallback angles are calculated, not manually specified.
Awarua's runner accepts polygon IDs but requires one shared simple source route;
the core currently picks the shortest gap between a supplied polygon pair.
Automatic network discovery, all-gap enumeration, headwater inference,
tributary junction handling and terrain/slope suitability checks are not yet
implemented in this experimental connector workflow.

## Jamestown connected-river heightmap trial

```sh
uv run python experiments/river_connectors/network_trial.py
uv run python experiments/river_connectors/build_heightmaps.py
uv run python experiments/river_connectors/review_network.py
```

Discovery now enumerates cached routes programmatically, merging degree-two LINZ
segments and splitting outside river/lake polygons. It requires at least one
**river** polygon anchor; lake-only tributaries and unattached networks are not
included. Different-river-polygon gaps use measured endpoint widths and linear
interpolation. Gaps that leave/re-enter the same polygon are excluded as possible
source drift. All original river water geometry remains authoritative.

For one-ended routes, the current interpretation checks the *free line end*:
if it touches a lake polygon or lies below sea level in the cached world DEM,
keep the receiving-river width constant. Otherwise, a free, unbranched endpoint
inside the playable footprint is treated as a headwater, starting at 0.5 m and
tapering linearly to the receiving river. This is a topology heuristic, not a
surveyed headwater or flow-direction claim. Mapped junctions and endpoints outside
the playable footprint are recorded for review, not assumed to be headwaters.
Only direct river-anchored routes are generated; unanchored branch trees await
a separate decision. Sea contact is DEM-mask based, not a paid bathymetry source.

Side tributary width sampling uses the mean of shortest complete cross-sections
inside the receiving polygon, searching orientations every 5 degrees. Using a
normal to the incoming tributary can otherwise measure lengthwise down the main
river. This is a local width estimate, not a surveyed tributary width or flow
split; broad receiving rivers still produce broad tributary mouths. Side entries
retain their incoming tangent and an interior anchor, rather than inferring a
direction from potentially coincident receiving-river cross-section midpoints.
Two-ended connectors retain the previously trialled sampling/alignment method.

Tributaries are calculated in full but their synthetic strips require width
**>=5 m** and absolute longitudinal terrain grade **<15 degrees**. Grade uses a
20 m centred chord along a raw-DEM profile smoothed over 10 m. Each failed reach
is omitted, including steep reaches encountered after an earlier passing reach.
The raster build independently applies those station-based thresholds at pixel
centres. No such extra cliff filter is applied to authoritative river polygons
or two-ended connectors in this pass. Smoothing is 20 m, addition-only, within
60 m of joins. Lake beds are excluded from new cuts.

Geometry output: `output/jamestown-river-lab/network-trial/` with generated and
skipped candidate reports, all widths/grades/stations, and geographic patches.
Heightmap output: `output/jamestown-connected-rivers/` containing paired 4096²
16-bit PNGs, a ZIP, map previews, geometry overlays, source provenance and audits.
The existing polygon-only package remains untouched. The build recomputes the
combined bed in new water plus a 60 m existing-river collar, removing old cutoff
banks without adding depths together. Lower absolute elevation targets win.
Synthetic patches are repaired/snapped to 1 mm after geographic round-trip to
avoid floating-point slivers; source polygons retain their coordinates.

Reference water levels remain estimates from the 30 m smoothed uncarved DEM;
the established lower Hollyford lake/sea reference is retained inside its source
polygon. Thus a 10 m synthetic depth cap is not a 10 m maximum change from the
previous terrain. Existing-valley detection and guaranteed downhill hydraulic
flow remain deferred. World-centre pixels are synchronised from the playable
map; no automatic in-game water sources are created.

### Current accepted export scope: reviewed connectors only

The broad inferred-tributary trial was rejected after review: Hollyford tributary
mouths were too wide. The preferred export now retains **only two Awarua and five
McKenzie polygon-to-polygon channels**. All newly inferred tributaries (including
McKenzie and Hokuri) and the added sea link are excluded. The existing polygon,
lake and seabed baseline is retained.

```sh
uv run python experiments/river_connectors/select_reviewed.py
uv run python experiments/river_connectors/build_heightmaps.py
uv run python experiments/river_connectors/review_network.py
```

`select_reviewed.py` filters the existing network trial with an explicit polygon
pair allowlist and verifies seven channels. It preserves the computed geometry
and writes `output/jamestown-river-lab/reviewed-connectors/`. The heightmap builder
now defaults to that reviewed selection; the earlier network geometry remains
available for research but is not implicitly exported. The regenerated
`jamestown-connected-rivers` package replaces the broader trial at the same output
path. Use `network_trial.py` then `select_reviewed.py` if regenerating from sources.
