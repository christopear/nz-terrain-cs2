# Jamestown river lab: lower Hollyford first

This branch isolates river work so we can inspect actual channel footprints and
original terrain before designing replacement excavation rules. Its initial
scope is the lower Hollyford between Lake McKerrow and the sea, within the current
Jamestown playable footprint. The map boundary splits it into two segments;
the lab does not connect them with an invented straight line.

## Current polygon-carver experiment

The preferred river profile is saved in
[`water-configs/river-preferred.json`](../water-configs/river-preferred.json):
10 m cap, 2 m shoreline transition, 5% initial gradient, 100% nominal interior
gradient and 1 m local crease smoothing. The tanh cap softens the actual slope.
Depth depends on distance to retained banks, not width normalisation. Existing
lake, lagoon and seabed configurations are unchanged.

The profile and smoothing now live in `water_carver.river`. The Jamestown
polygon-only trial uses them over a matching cached sea/lake baseline:

```sh
uv run python -m nz_terrain.river.polygon_build --baseline-root /path/to/baseline-checkout
uv run python -m nz_terrain.river.review_package --baseline-root /path/to/baseline-checkout
```

Output: `output/jamestown-polygon-rivers/`, with paired heightmaps, previews,
metadata and a ZIP. The ordinary `--with-rivers` command still uses the legacy
centreline pipeline; use the command above for this trial. It verifies the
baseline location and water settings, including presence of Waiuna Lagoon.

The playable and world grids use fourfold subpixel depth sampling, then average
water-covered samples. Only polygon-centre pixels are applied; land outside the
river footprints is unchanged. Both grids retain their existing encoding and
rotation. Changed playable cells are averaged into the world map centre.
Hollyford uses a synthetic reference from 1 m at the lake to 0 m at the sea;
other polygons use 30 m masked smoothing of their original DEM surface. Absolute
target beds are combined by taking the lower elevation, never adding depths.
The two reviewed Hollyford connections feather into existing water over 100 m,
with a 1.5-pixel collar at each opening to close DEM/raster classification gaps.

Other polygon end caps, river-stage estimation, steep reaches, and existing
valley detection still need review. This is an in-game trial, not a finished
hydrology model. The historical centreline lab below is reference material.

### Hollyford openings

Reproduce the diagnostic using a checkout containing cached lake metadata:

```sh
uv run --all-packages python packages/water-carver/examples/inspect_hollyford_openings.py --baseline-root /path/to/baseline-checkout
```

Output: `output/jamestown-river-lab/openings/`, containing a before/after plot,
opening candidates, retained bank geometry and an evidence report. The script
reads the saved preference and leaves production inputs untouched.

Two exterior route crossings identify the expected cutoffs:

- The 121.1 m lake-end edge exactly coincides with Lake McKerrow's boundary.
- The 174.4 m mouth edge crosses the route at 87.5 degrees, 77.2 m before its
  downstream endpoint. This remains a candidate until checked against a sea mask.

Exclude these edges from bank-distance calculation, while retaining the polygon
as the water footprint. All island boundaries remain banks. The diagnostic
checks both retained island lengths and the total excluded boundary length.
The route also crosses an island twice: removing every route crossing would
therefore be unsafe. Straight edges alone are not evidence of openings either.

For broader automation, use shared boundaries with adjoining water polygons as
primary evidence, with exterior route crossings, continuation and crossing angle
as supporting evidence. Review ambiguous candidates. Derive banks before map
clipping so map edges do not become dams. Different datasets may need a careful
geometry-matching tolerance; this experiment uses exact shared boundaries.

The opening diagnostic plots assume a flat water surface and only remove
artificial river-end banks. The separate map build above adds the lake/sea
connections and estimated water elevations for in-game review. A broad flat bed in wide reaches is
still an accepted limitation of the selected bank-distance profile.

## Historical centreline lab

## What exists in this first iteration

- `src/nz_terrain/river/core.py`: existing profile solver, terrain sampler and
  segment carver extracted from the map-package orchestration. These are legacy
  algorithms, unchanged; extracting them does **not** solve the reported issues.
- `rivers.py`: retains its existing public functions and commands through a small
  adapter. Existing map generation continues to use the legacy behaviour.
- `river/geometry.py`: independent metre-based transects and polygon widths.
- `river/lab.py`: a Jamestown-only diagnostic command using the **original,
  uncarved** cached land surface, existing Hollyford route and mapped river-area
  polygon. Existing lake boundaries exclude lake interiors; no lake or ocean
  builder is invoked.

The lab produces CSV profiles, review plots, a GeoJSON transect layer and a JSON
report. It does not carve terrain, regenerate packages, or place water sources.
In particular, variable-width excavation, steep-section exclusion, source basins
and existing-valley preservation are **not implemented yet**. The 45-degree
threshold currently produces an along-route diagnostic flag only, evaluated at
the configured station spacing; short cliffs between stations can be missed.

## Running it

```sh
uv sync --locked

# With an existing matching Jamestown baseline/cache in this checkout:
uv run nz-terrain-river-lab

# Read an existing checkout's baseline/cache from a separate worktree:
uv run nz-terrain-river-lab --baseline-root /path/to/existing/nz-terrain

# Optional: explicitly replace the committed polygon snapshot from LINZ:
uv run nz-terrain-river-lab fetch-polygons

# Numerical/geometry regression tests, no downloads needed:
uv run python -m unittest discover -s tests -v
```

The default diagnostic destination is `output/jamestown-river-lab/` in the current
checkout. `--baseline-root` is a read-only input; outputs stay in this checkout
unless `--output` is supplied. No raster download is needed for inspection once
the existing baseline is available. A location/grid mismatch fails explicitly.

- `stations.csv`: varying mapped spans, original centreline elevations,
  along-route slope and unknown-width reasons.
- `terrain-cross-sections.csv`: original DEM elevations at 5 m cross-section
  intervals; these are not underwater bed surveys.
- `transects.geojson`: section locations for GIS review, with source attribution.
- `width-review.png` and `terrain-cross-sections.png`: visual review.
- `report.json`: settings, source hashes, width summaries and limitations.

Widths use a local metre projection, perpendicular to a 50 m route-tangent
window. Only the connected polygon interval containing the route point counts.
An island splits that interval; disconnected braids are not summed. A centreline
outside the polygon, ambiguous intersection, or a transect too short to reach
both banks yields **unknown**, never zero or a default 5 m channel. Off-grid DEM
samples remain unknown rather than repeating the nearest edge elevation.

The supplied 27 October 2026 snapshot and default 50 m station spacing produced
121 stations: **117 mapped spans, 4 outside-polygon cases**. Mapped spans were
about **55–403 m, median 179 m**. This is a diagnostic measurement of the polygon,
not a field survey or width at a known water stage. Bends, islands, bank
alignment and mouth geometry require review before converting widths to a bed.

## What the source geometries contain

LINZ distinguishes [river_cl](https://docs.topo.linz.govt.nz/data-dictionary/tdd-class-river_cl.html)
from [river_poly](https://docs.topo.linz.govt.nz/data-dictionary/tdd-class-river_poly.html).
The first is a line representation for rivers too small to map as areas at
1:50,000; the second is an area representation for wider rivers. A line traces a
route; a polygon represents a footprint, potentially with holes and multiple parts.

On 27 October 2026, both public pilot services reported `hasZ=false` and
`hasM=false`: their geometry is 2D. Neither schema had a measured width, bed-depth,
water-level or discharge field. Polygon widths must be derived from geometry;
line widths require other evidence. Other GIS products can contain Z coordinates
or width attributes, but Z alone does not establish whether it represents the
land, water surface or river bed.

- [River name lines service](https://services.arcgis.com/xdsHIIxuCWByZiCB/ArcGIS/rest/services/LINZ_NZ_River_Name_Lines_Pilot/FeatureServer/0)
- [River name polygons service](https://services.arcgis.com/xdsHIIxuCWByZiCB/ArcGIS/rest/services/LINZ_NZ_River_Name_Polygons_Pilot/FeatureServer/0)

The committed `river-configs/jamestown-hollyford-polygons.geojson` retains the
service URL, query URLs, retrieval time, location, service schema and attribution.
It is LINZ data under CC BY 4.0; the code's MIT licence does not replace that licence.
The existing Hollyford route is OSM-derived and retains its ODbL attribution.

## Width evidence worth pursuing

1. **Mapped river polygons first.** They already describe broad, varying channel
   footprints for this river. Keep their outline and islands, review alignment,
   then define how much of that footprint should be a flat bottom versus banks.
2. **Free georeferenced aerial imagery for verification and gaps.** LINZ publishes
   [free aerial imagery](https://www.linz.govt.nz/products-services/data/types-linz-data/aerial-imagery),
   commonly much finer than satellite imagery. Trace wet banks manually for the
   first test reach; later consider assisted segmentation with reviewed masks.
   Clouds, shadows, vegetation, shallow gravel and capture-stage changes make
   automatic classification uncertain. Aerial imagery cannot by itself supply
   reliable riverbed depth.
3. **Satellite water masks for broad channels and time variation.**
   [Sentinel-2](https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-2/Facts_and_figures)
   has 10 m visible/NIR bands. That is useful scale for a roughly 200 m river,
   but a 5–10 m stream is about one pixel or less: not reliable bank geometry.
   Multi-date water masks could distinguish persistent water from occasional
   inundation. Use original multispectral imagery for water-index analysis, not
   an RGB-only visual basemap.
4. **Existing width datasets and local survey reports as checks.**
   [GRWL](https://water.usgs.gov/catalog/datasets/120270a9-e0b6-42d8-9b1f-17db852fd2b4/)
   derives river widths from 30 m Landsat imagery. It is a useful broad-river
   cross-check, not a precision source for small streams. Hollyford coverage
   has not been verified. Council cross-sections or reports could provide local
   widths, but no such Hollyford survey has been acquired for this PR.
5. **Manual width anchors where evidence is missing.** Record source, date and
   uncertainty, then interpolate between reviewed anchors. These are explicit
   design decisions, preferable to presenting stream order or upstream length
   as a measured width.

### Confirmed Hollyford aerial coverage

The public LINZ imagery STAC index was checked, not just searched by region name.
Two tiles from the **Southland 0.75 m 2005–2011** collection cover the complete
Hollyford polygon used here:

- [CA08_10000_0105](https://nz-imagery.s3.ap-southeast-2.amazonaws.com/southland/southland_2005-2011_0.75m/rgb/2193/CA08_10000_0105.json)
- [CA08_10000_0205](https://nz-imagery.s3.ap-southeast-2.amazonaws.com/southland/southland_2005-2011_0.75m/rgb/2193/CA08_10000_0205.json)

Both items expose the collection-era capture interval, not an exact capture day.
The 2026 metadata processing date is **not** the imagery capture date. The modern
[Southland 2023–2024 0.25 m collection](https://nz-imagery.s3.ap-southeast-2.amazonaws.com/southland/southland_2023_0.25m/rgb/2193/collection.json)
does not extend north to this river; the checked West Coast 2024–2026 collection
ends north/east of it. This is not an exhaustive search for every possible newer
image. Older aerial imagery should be compared with newer satellite observations
before treating a mobile mouth or bank as current.

## Next iteration to agree together

Choose whether the target width is ordinary wetted channel, bankfull channel, or
a deliberately simplified game channel. These are different quantities. Then:

1. Review polygon/route alignment and width outliers on this single reach.
2. Inspect cross-sections to identify existing usable channel space. Design the
   minimum needed excavation, not a fixed subtraction from a valley floor.
3. Define tapering and bank slopes using reviewed geometry/anchors.
4. Define slope exclusion over a stated metre window, restart behaviour below
   steep sections, and source-basin constraints. Keep cliffs uncarved; do not
   create a basin where doing so would undermine a bank or require excessive cuts.
5. Emit an independent absolute river-bed target plus a modification mask and
   decision report; compose with fixed sea/lake boundary conditions only later.

No other locations should be migrated to the replacement carver until the
Jamestown test is reviewed in-game.


## Additional polygon-trial presets

Okarito, Jackson Bay and Takaka use the same independent river profile and
current lagoon defaults. Each has a saved LINZ river-area snapshot and explicit
reviewed opening edges under `river-configs/<location>-{polygons,openings}.geojson`.
Snapshots retain full polygons returned for the playable geographic envelope;
they do not claim complete coverage of the surrounding world.

```sh
uv run python -m nz_terrain.river.fetch_polygons --location okarito
uv run nz-terrain --location okarito
uv run python -m nz_terrain.river.polygon_build --location okarito --baseline-root .
uv run python -m nz_terrain.river.review_package --location okarito --baseline-root .
```

Substitute `jackson-bay` or `takaka` as appropriate. Downloads are a separate
explicit step; a fetch that changes polygon geometry requires reviewing existing
opening annotations before building. Named presets, current baseline settings,
and the original uncarved surface cache must match.

Openings in these builds include exact shared lake interfaces (within a 5 cm
matching tolerance, with an outward probe inside the lake), the 120 m Takaka
estuary cutoff and the 325 m Arawhata mouth cutoff. Other ambiguous crossings
remain banks. Reference river levels use the same 30 m smoothed original DEM
estimate as the Jamestown trial outside Hollyford. They are not a hydraulic
surface solver. No centreline fallback is applied; notably the downloaded
Takaka river-area layer has negligible coverage along the Waikoropupu route.

Bank distances use an exact spatial index of original bank segments for larger
polygons. A regression test compares indexed distances, the resulting profile
and smoothing with direct vector distances, including an island and an opening.
The full source geometry is retained; there is no polygon simplification.
