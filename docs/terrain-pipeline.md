# nz-terrain

Generate separate, true-scale vanilla CS:II terrain packages from named locations.

```sh
uv sync
uv run nz-terrain list
uv run nz-terrain --location jamestown
uv run nz-terrain --location takaka --location port-underwood
uv run nz-terrain --all
```

For footprint polygons without downloading elevation data:

```sh
uv run nz-terrain footprints --all
```

## Locations

Edit `locations.json` to add or change a preset. Coordinates are **longitude, latitude**.

| Slug | Longitude | Latitude | Map-up bearing |
|---|---:|---:|---:|
| jamestown | 168.09214 | -44.35322 | 345° |
| takaka | 172.80447 | -40.85758 | 30° |
| port-underwood | 174.05366 | -41.30641 | 30° |
| bay-of-islands | 174.09303 | -35.26098 | 0° |
| milford-sound | 167.95738 | -44.69578 | 0° |
| jackson-bay | 168.66856 | -44.01286 | 0° |
| okarito | 170.22337 | -43.21401 | 305° |

**Rotation means the direction of the top of the map, clockwise from true north.** Takaka's top edge faces bearing 030°. The sampling footprint is rotated in metre coordinates; the exported image is never rotated/resized afterward. Every map retains the requested centre and true-scale dimensions. A local transverse Mercator projection has unit scale at the centre with small distortion away from it.

The old Jamestown centre is retired.

## Output and caching

Each completed location has its own folder:

```text
output/takaka/
  takaka-heightmaps.zip
  IMPORT.txt
  heightmaps/
    takaka-playable.png
    takaka-world.png
  previews/
    playable-preview.jpg
    world-preview.jpg
  metadata/
    footprints.geojson
    manifest.json
```

The ZIP preserves this same folder layout. Copy only the two PNGs from `heightmaps/` into the game’s Heightmaps folder. `previews/` contains reference images; `metadata/` contains manifests, geographic data and, for river candidates, reports and profile CSVs. `IMPORT.txt` stays at the top level alongside the ZIP. PNG names include the location slug so several locations can coexist in the game’s Heightmaps folder. Packages are published only after generation and validation succeed. A failed build leaves the previous completed package intact.

Source grids live under `data/locations/<slug>/<config-hash>/`; changes to centre, rotation or source lists use a new cache. Shared LINZ metadata lives in `data/catalogue/`. Completed source grids are reused on reruns. To refresh source metadata, remove `data/catalogue/` and the relevant location cache explicitly. Cached metadata is a reproducible snapshot, not an automatically refreshed catalogue.

`--workers 3` is the default for concurrent tile downloads. Lower it if memory/network capacity is limited. Each in-flight grid is about 64 MB; a build may use roughly 1–2 GB RAM during blending and seabed processing. Network access is required for uncached source windows. Local output/source data is ignored by Git.

## Sources

The pipeline selects actual intersecting LINZ STAC items. It uses Topo50 tile names as a fast spatial prefilter, then checks STAC bounding boxes before reading rasters. Recognised tile-name geometry is covered by tests; unknown conventions are fetched for direct inspection.

Source priority:

1. Configured **LINZ bare-earth LiDAR DEM** collections, oldest to newest. Later entries override overlapping older entries. Configure their paths in `lidar_collections`; the program does not promise to discover every newly published survey automatically.
2. **LINZ NZ 8 m contour DEM**, where covered. Its cells are interpolated mainly from Topo50 20 m contours; this is not 8 m surveyed detail.
3. **Copernicus GLO-30 Public** where LINZ is absent. It represents surfaces including vegetation.

LiDAR uses published 2 m overviews for the 3.5 m playable grid and 8 m overviews for the 14 m scenery grid where available. Outer LiDAR coverage is feathered over 250 m. Source elevations are not independently normalised. Cross-source vertical datums are **not formally harmonised**: this is a cartographic/game map, not a survey product.

- [LINZ public elevation catalogue](https://nz-elevation.s3.ap-southeast-2.amazonaws.com/catalog.json)
- [LINZ elevation archive documentation](https://github.com/linz/elevation)
- [LINZ 8 m DEM metadata](https://data.linz.govt.nz/layer/51768-nz-8m-digital-elevation-model-2012/)
- [Copernicus public dataset](https://registry.opendata.aws/copernicus-dem/)

Ocean-only Copernicus tiles are omitted by the publisher. A 404 is treated as
ocean only after verifying the tile is absent from both GLO-30 and the worldwide
GLO-90 catalogues; the fill is confined to that tile and recorded in the source
manifest. Listed-but-unavailable tiles and other missing coverage still fail the
build. This handles offshore scenery corners at Okarito without
discarding land.

Collections and individual source URLs are retained in each manifest. Attribution is included in each package. Test fixtures are small excerpts from the public LINZ STAC metadata, licensed CC BY 4.0.

## Synthetic sea, lake and lagoon beds

All underwater geometry is synthetic. No paid bathymetry is used. The pipeline
keeps the original blended land surface in the cache, calculates independent
**absolute bed elevations**, then combines their domains. Depths are never added
or subtracted repeatedly. Where targets overlap, the lower bed wins; averaging
would risk raising a riverbed into a sill at a lake entrance.

```sh
uv run nz-terrain --location jamestown --with-rivers \
  --default-lake-depth=30 --default-lake-gradient=0.05 \
  --default-lagoon-depth=10 --default-lagoon-gradient=0.05
uv run nz-terrain --all --with-rivers
```

`--with-rivers` also rebuilds the separate river candidate wherever a river
configuration exists (all configured presets). Takaka retains its four manually tuned
main rivers; the other presets use mapped-network selection. Rebuild river candidates after changing water settings
so both stages refer to the same baseline.

| Domain | Default depth ceiling | Near-shore gradient | Water-surface reference |
|---|---:|---:|---|
| Sea | 30 m | 0.05 (5%) | 0 m |
| Lake | 30 m | 0.05 (5%) | Per-lake elevation |
| Lagoon | 10 m | 1% shoreline ramping to 5% over 25 m | Per-lagoon elevation |

Sea and ordinary lakes use `depth = D × tanh(distance × gradient / D)`.
Lagoons use a smooth gradient transition from 1% at the shore to a nominal 5%
over 25 horizontal metres, then the same tanh depth cap. Distance smoothing is
25 m for lagoons and remains 14 m for lakes and sea. These are bed slopes, not
water-surface gradients or flow settings. Island shorelines are retained; narrow
water bodies may be much shallower than the ceiling. The default sea formula
remains equivalent to `30 × tanh(distance / 600)`.

Lagoon defaults are `lagoon_distance_smoothing_m: 25`,
`lagoon_shore_transition_m: 25` and `lagoon_shore_gradient: 0.01`. A zero
transition restores the plain tanh profile. Per-feature overrides can specify
`distance_smoothing_m`, `shore_transition_m` and `shore_gradient`.

All six `--default-{sea,lake,lagoon}-{depth,gradient}` options are available.
Persistent defaults are in `water-configs/defaults.json`; optional location files
such as `water-configs/jamestown.json` can override them. CLI values override
location defaults. Per-feature overrides, keyed by LINZ `t50_fid`, take priority
and can specify `depth_m`, `dropoff_m`, `water_level_m`, `kind`, or `skip`.
A per-feature `dropoff_m` fixes the horizontal scale rather than deriving it from
the default gradient. This keeps McKerrow's sourced depth parameter and chosen
600 m transition separate from generic defaults.

Lake/lagoon polygons explicitly exclude those interiors from the sea pass,
including lakebeds below sea level. Polygon holes preserve islands. Lakes are
synthesised on both grids with padded distance calculations to avoid false shores
at map edges; affected world-centre cells use 4×4 playable averages. The original
land surface is preserved outside the relevant water domains.

Boundaries come from free LINZ Topo50 lake polygons and a public LINZ-derived
lagoon layer. Names containing “Lagoon” also receive lagoon settings—Waiuna is
classified as a lake in the source dataset. Seasonal/dry lakes are excluded.
Missing/unreliable surface-elevation attributes fall back to the median original
DEM elevation and are flagged as estimates. Coverage is incomplete: no returned
features does not prove no lakes, where source coverage is incomplete. Separate
pond layers and unmapped water are not automatically included.

McKerrow uses a published approximate maximum-depth parameter of **109 m** and
LINZ's **+1 m** surface level. Its basin shape remains entirely synthetic; the
achieved depth depends on the footprint and distance from shore. Waiuna uses the default 10 m lagoon ceiling, 5% main gradient, 25 m shoreline
transition and 25 m distance smoothing. With the source's +15 m water level,
its generated maximum depth is about 9.24 m (bed elevation about +5.76 m).
See each package's `metadata/water-report.json` for actual parameters, source
links, inferred levels, exclusions and per-grid depths.

- [LINZ lake definitions and water-level attributes](https://docs.topo.linz.govt.nz/data-dictionary/tdd-class-lake_poly.html)
- [Lake McKerrow depth reference, Wellman and Wilson (1965)](https://www.tandfonline.com/doi/pdf/10.1080/00288306.1965.10426432)

These packages contain terrain, not water sources. The editor still needs water
at the reported lake/lagoon levels and suitable river inflows. Sources, outlet
behaviour and flooding have not been calibrated in CS2.

## Mapped river networks

```sh
uv run nz-terrain-rivers --location jamestown --minimum-river-width=5 \
  --default-river-depth=1.5
```

Network configurations live in `river-configs/<slug>.json`. Prepare mainland
line snapshots explicitly with `uv run nz-terrain-discover-rivers --all`; this
preserves existing OSM supplements. Wide channels omitted by the LINZ line layer
are supplemented from OSM/WaterwayMap; source URLs are stored on those features.
Rebuild baseline and rivers together after refreshing water data.

Takaka keeps its existing four tuned channels. Other maps use the 5 m estimated
width threshold below, including trimming each river's narrow upstream sections.
See each candidate's report for included reaches and source/terrain exclusions.

The Jamestown network configuration is `river-configs/jamestown.json`. The first broad pass
retains over 100 reaches, including Awarua, Dry Awarua, the in-footprint lower
Hollyford, McKenzie/Jerusalem/Hokuri/Wilmot sections, and unnamed streams around
Jamestown/Gravel Cove. The separate Kaipo River is west of the playable footprint;
**Hollyford is McKerrow's outlet**. River geometry is clipped to the playable
footprint; portions outside it are not excavated in the scenery map.

The **5 m selection threshold trims narrow upstream portions, not entire rivers**.
Widths are not recorded on these centreline datasets. The current explicit
heuristic is `5 × sqrt(upstream mapped channel kilometres)`, with assumptions at
boundary/lake inflows and named overrides. It is an experimental size proxy,
not a measured five-metre cutoff. Retained beds use at least 7 m bottoms (two
heightmap pixels) for reliable rasterisation. Adjust `width_coefficient` and
`minimum_river_width_m` independently of the carved width/depth parameters.

Connected endpoints share bed elevations; bounded snaps are recorded. Mapped
lake/lagoon interiors are excluded from river routing, and adjoining river
surface proxies blend to the water body's level over 200 m. Trusted OSM flow
direction takes precedence over DEM endpoint height; other components infer
outlets from their lowest endpoint. Incomplete source geometry can leave gaps.
The lower Hollyford has an explicit +1 m surface proxy, 3 m bed depth and a local
20 m excavation cap because the coarse land DEM places false raised ground along
its outlet. Other network reaches default to a 10 m cut cap. Preferred downhill
slope is soft; flat lake outlets are allowed.

Reaches that cannot meet the excavation limit retain their longest feasible
downstream portion; excluded upper lengths are reported. Cyclic side channels,
if present, are reduced to a spanning tree and counted. Check
`metadata/river-report.json` for every retained/excluded reach, width-based trim,
profile-based trim, source marker and raster-level uphill-step warning.
Source markers are placement candidates, not instructions to add an independent
inflow to every confluence or lake outlet. Filling every minor stream requires
in-game tuning and may be unnecessary for the landscape.

## Game format and import

- Playable terrain: **14,336 m × 14,336 m**, 3.5 m per pixel.
- Scenery: **57,344 m × 57,344 m**, 14 m per pixel.
- Both: **4096 × 4096 unsigned 16-bit grayscale PNG**.
- Encoding: `elevation_m = pixel / 65535 × 4096 − 512`.
- No vertical exaggeration. Black is −512 m, white is 3584 m.

Use the same settings as the successful Jamestown import. A zero-based 4096 m import span places real sea level at game elevation 512 m. Check your editor's actual scaling/offset controls; the generated packages have not been tested in-game by this project. Use both PNGs from the same location/package.

Files normally go in `%USERPROFILE%\AppData\LocalLow\Colossal Order\Cities Skylines II\Heightmaps`. See each package's `IMPORT.txt`. Importing requires no mods. Do not auto-contrast, resize or convert these images to 8-bit.

## Experimental Takaka river beds

```sh
uv run nz-terrain --location takaka
uv run nz-terrain-rivers --location takaka
```

The second command uses the completed baseline, checks its location matches,
and writes a separate `output/takaka-rivers/takaka-rivers-heightmaps.zip`.
Use both `takaka-rivers-playable.png` and `takaka-rivers-world.png` together.
The ordinary Takaka package remains available for comparison. No mods required;
water sources still need manual placement in the editor.

Edit `river-configs/takaka.json` for independent bottom width, bank-blending width,
excavation depth, preferred slope and maximum cut. These are experimental design
choices, not surveyed dimensions or guaranteed simulated water depths:

| River | Carved length | Bottom width | Full edit width | Target excavation below surface proxy |
|---|---:|---:|---:|---:|
| Takaka | 15.82 km | 18 m | 48 m | 3 m |
| Waingaro | 6.21 km | 14 m | 36 m | 2.5 m |
| Anatoki | 10.85 km | 12 m | 32 m | 2.5 m |
| Lower Waikoropupu | 5.71 km | 10 m | 26 m | 2 m |

The route database in `river-configs/takaka-routes.geojson` is editable WGS84
geometry derived from OpenStreetMap via WaterwayMap. Attribution and download
links are in `river-configs/SOURCES.md`. The disconnected mountainous Anatoki
source fragment is excluded; all selected lower reaches connect to Takaka.
The first 4 km of the clipped mountainous Waikoropupu route are deliberately
unchanged: a >12 m uphill excursion in the sampled surface made the initial
profile infeasible within the excavation limit. Its lower valley is included.
This is a first experiment, not complete catchment conditioning. Springs River
and unconfigured branches are unchanged; synthetic lakes/lagoons are inherited
from the baseline package.

The algorithm samples a small corridor around each line, smooths a surface proxy,
then solves the longitudinal bed with hard downhill and cut-depth constraints.
Tributary outlets share their main-river junction elevation. A preferred 0.25%
bed slope is **soft**: a fixed sea outlet and limited natural drop take priority.
The main reach has only about 35 m of surface fall over 15.8 km; a steeper uniform
grade cannot be manufactured while preserving both endpoints. The outlet bed is
−3 m, which does not change the sea's water level. Steep natural reaches can retain
their fall; there is no automatic waterfall detection or forced maximum gradient.

Each profile targets a flat central bed and cosine-blended banks using the
original source surface. It combines with baseline water beds using the lower
absolute elevation, so it cannot add its depth to a previously excavated lake. Every pixel has a 10 m excavation cap. Outside these corridors the
playable PNG is bit-for-bit unchanged. The world map's affected central cells are
updated from 4×4 playable averages; river carving ends at the playable boundary.
Sources therefore need to be inside the playable area, not in the scenery map.

The ZIP includes route overview, longitudinal profiles and example cross-sections
in `previews/`, plus profile CSVs and `river-report.json` in `metadata/` with source markers and numerical checks.
Start with modest inflow at the markers, let the simulation settle, and adjust
until the narrow channel holds roughly 1.5–3 m of water without overtopping.
Inspect all three junctions and the estuary. Change source supply and width
before imposing more bed slope. Source strength is not calibrated and no real
flow rate is being equated to a CS2 editor setting. Terrain excavation volume is
reported in cubic metres; it is neither water volume nor discharge.

Council flood modelling also required surveyed river bathymetry beyond the
LiDAR: [Takaka River Modelling Report, April 2024](https://www.tasman.govt.nz/document/serve/Takaka%20River%20Modelling%20Report%20April%202024.pdf?DocID=36203).
This project has not obtained that survey or run an in-game/hydraulic simulation.

## Validation

```sh
uv run python -m unittest discover -s tests -v
uv run nz-terrain-validate
```

Tests check rotation/scale, source tile bounds, depth profiles, encoding limits, island preservation, absolute-target composition/idempotence, partial headwater trimming and trusted flow direction despite misleading DEM elevations. Each build also checks complete fallback coverage, rejects clipping, reopens both PNGs to verify uint16 dimensions/pixels, asserts unchanged land/shoreline during seabed rebuilding, checks fully submerged centre alignment to within half an encoding unit (~3.13 cm), and verifies the ZIP contents. Mixed shoreline cells prioritise land preservation over exact averaging.

`nz-terrain-validate` checks all configured baseline/candidate pairs against the current
location presets, verifies checksums and archive contents, and confirms river
carving only lowers the playable terrain and changes the corresponding scenery
cells. It writes `output/validation-summary.json` and the contact sheet
`output/locations-preview.jpg`. River obstructions remain explicit report flags;
passing file checks does not imply successful in-game water simulation.
