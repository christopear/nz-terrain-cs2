# nz-terrain-cs2

Generate true-scale New Zealand terrain packages for **Cities: Skylines II**, using
free public elevation data and synthetic sea, lake, lagoon and river beds.
No mods or paid bathymetry are required. Water sources must be placed in-game.

## Requirements and setup

- Python **3.11 or newer**.
- [uv](https://docs.astral.sh/uv/getting-started/installation/) to install Python
  dependencies and run commands.
- Internet access for uncached public GIS data, several GB of disk space, and
  roughly 1–2 GB of RAM per terrain build (more with higher concurrency).

From the repository directory:

```sh
uv sync --locked
uv run nz-terrain list
uv run nz-terrain --location jamestown --with-rivers
```

`uv sync` installs this project as an editable package and uses `uv.lock` for
reproducible dependency versions. Native GIS dependencies normally come from
Python wheels on supported platforms; platforms without compatible wheels may
need system GDAL/PROJ libraries and a compiler.

## Usage

```sh
# Build one baseline terrain package, or several presets
uv run nz-terrain --location milford-sound
uv run nz-terrain --location takaka --location port-underwood

# Build all presets and their configured river candidates
uv run nz-terrain --all --with-rivers

# Export footprints without downloading elevations
uv run nz-terrain footprints --all

# Set synthetic water-bed defaults (gradient is metres per metre)
uv run nz-terrain --location jamestown --with-rivers \
  --default-lake-depth 30 --default-lake-gradient 0.05 \
  --default-lagoon-depth 10 --default-lagoon-gradient 0.02

# Rebuild rivers against an existing, matching baseline
uv run nz-terrain-rivers --location jamestown --minimum-river-width 5
```

`uv run python -m nz_terrain` is equivalent to `uv run nz-terrain`.
Use `--help` on the terrain and river commands for options. Run commands from
this checkout. If using an installed wheel elsewhere, set `NZ_TERRAIN_ROOT` to
the checkout containing the presets and data directories; configuration files
and GIS snapshots are maintained in the checkout, not bundled into the wheel.

The seven presets are **Jamestown, Takaka, Port Underwood, Bay of Islands,
Milford Sound, Jackson Bay and Okarito**. Their centres,
clockwise map-up bearings and elevation collections are in `locations.json`.
See [the pipeline guide](docs/terrain-pipeline.md) for the coordinate table,
water parameters, river selection rules and import details.

## Outputs

Each map has its own package. A river candidate is generated separately from
its baseline:

```text
output/jamestown/                # Synthetic sea, lake and lagoon beds
output/jamestown-rivers/         # Same baseline plus selected river beds
  jamestown-rivers-heightmaps.zip
  IMPORT.txt
  heightmaps/                    # Playable and scenery PNGs to import
  previews/                      # Terrain, river routes and profiles
  metadata/                      # Source manifests, checksums and diagnostics
```

Both PNGs are **4096 × 4096 unsigned 16-bit grayscale**. The playable map is
14,336 m square; scenery is 57,344 m square. Use both PNGs from the same package.
The elevation encoding is `metres = pixel / 65535 × 4096 − 512` with no vertical
exaggeration. Follow `IMPORT.txt`; do not resize or auto-contrast the PNGs.

Downloaded metadata and source grids are cached under `data/`. Changes to preset
coordinates, rotation or source collections select a new grid cache. Cached
source metadata is a snapshot, not automatically refreshed live data. Generated
maps and caches are ignored by Git; the curated GeoJSON inputs are tracked.

## Data sources and limitations

Elevation priority is configured **LINZ bare-earth LiDAR**, then the **LINZ 8 m
contour DEM**, then **Copernicus GLO-30**. LiDAR collections currently require
explicit configuration. The 8 m product is interpolated terrain, and Copernicus
includes vegetation and structures. Vertical datums are not formally harmonised.

Lake/lagoon boundaries and river centrelines come from LINZ-derived public GIS
services, with OpenStreetMap supplements. All underwater shapes are synthetic.
A sourced approximate depth can constrain a lake; it does not supply surveyed
bathymetry. Sea, lake and river targets are calculated independently and combined
as absolute elevations, avoiding repeated subtraction of depths.

Read [data sources, attribution and acquisition](docs/data-sources.md) before
refreshing or redistributing data. Per-package manifests record actual sources.
Third-party data retains its own licence; this repository does not relicense it.

River widths are heuristic estimates. The 5 m cutoff trims narrow upstream
sections; retained beds have a minimum raster-compatible width. Excavation
limits may trim additional reaches. Reports identify exclusions and uphill
raster steps. Water flow and source strengths have **not been validated in CS2**.

### Can I use arbitrary NZ coordinates?

The underlying projection and raster sampling support other coordinates, but
there is currently **no arbitrary latitude/longitude CLI** or complete automatic
location setup. Add a preset to `locations.json`, choose suitable LiDAR
collections, prepare river geometry/configuration, and review water coverage and
profiles. Coverage gaps, coastal classification, lake levels and difficult
river terrain can still require local judgement. The current presets are
configured examples, not a guarantee of unattended results anywhere in NZ.

## Repository layout

```text
src/nz_terrain/    # Installable Python package and CLI implementations
locations.json    # Location centres, rotations and elevation collections
river-configs/    # Curated river routes, parameters and source attribution
water-configs/    # Synthetic basin defaults, overrides and extra polygons
tests/            # Unit tests and small source-metadata fixtures
docs/             # Data-source and detailed pipeline documentation
data/             # Ignored download/source-grid cache
output/           # Ignored generated packages and previews
pyproject.toml    # Project metadata, dependencies and CLI entry points
uv.lock           # Tracked application dependency lockfile
```

## Development and validation

```sh
uv sync --locked
uv run python -m unittest discover -s tests -v
uv run nz-terrain-validate
uv build
```

Unit tests require no downloads. The output validator requires completed baseline
and river packages for every configured preset. It checks their locations,
PNG dimensions/checksums, archive integrity and playable/scenery consistency,
then writes `output/validation-summary.json` and `output/locations-preview.jpg`.
File validation does not establish successful in-game water flow.

The `.gitignore` starts from [GitHub's Python template](https://github.com/github/gitignore/blob/main/Python.gitignore),
with project-specific cache, output and local-tool exclusions. Keep `uv.lock`
tracked when changing dependencies.

## Jamestown river development

The [lower-Hollyford river lab](docs/river-lab.md) isolates river diagnostics from
sea/lake generation. Run `uv run nz-terrain-river-lab` against an existing
Jamestown baseline, or pass `--baseline-root /path/to/existing/checkout` from a
worktree. It measures polygon-derived widths and original terrain cross-sections;
it does not yet implement the replacement carver or change map packages.
