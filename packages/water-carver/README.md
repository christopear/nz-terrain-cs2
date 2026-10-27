# water-carver

An independent Python/uv package for generating synthetic water beds from 2D
polygons. It has no dependency on `nz-terrain`, its datasets, or game encodings.
The existing NZ terrain pipeline does not use it yet.

## Setup and tests

From the repository root:

```sh
uv sync --all-packages --locked
uv run --package water-carver python -m unittest discover -s packages/water-carver/tests -v
uv build --package water-carver
```

Alternatively, the package directory can be used as a standalone project outside
this workspace with `uv sync`. A workspace shares the root `uv.lock`.

## Example

```python
from affine import Affine
from shapely.geometry import box
from water_carver import BedParameters, Grid, generate_bed

# All polygon and affine coordinates are already in NZTM metres.
water = box(1_400_000, 5_000_000, 1_401_000, 5_001_000)
grid = Grid(
    shape=(512, 512),
    transform=Affine.translation(1_400_000, 5_001_024) * Affine.scale(2, -2),
    crs="EPSG:2193",
)
bed = generate_bed(
    water,
    grid,
    water_level_m=20.0,
    parameters=BedParameters(max_depth_m=10, bank_gradient=0.05, profile="tanh"),
)
# bed.mask: cells inside water, sampled at pixel centres
# bed.depth_m: positive depth below the supplied reference water surface
# bed.elevation_m: absolute elevation = water_level_m - depth_m
# Both numeric arrays contain NaN outside water. No terrain is modified.
```

`water_level_m` can also be an array matching the grid, enabling a future caller
to supply a sloping river surface. It must be finite inside the water polygon.
The caller is responsible for this surface, its vertical datum and consistency
with any terrain to which the result is eventually applied. The carver does not
infer water levels from banks or existing depressions.

## Parameters and profiles

Let `d` be the shortest distance to a bank, `D` the maximum depth and `g` the
initial bank gradient in m/m (0.05 means 5%, **not** 5 degrees).

| Profile | Depth | Behaviour |
|---|---|---|
| `linear` | `min(g*d, D)` | Straight bank slope, then a flat depth cap |
| `tanh` (default) | `D*tanh(g*d/D)` | Smooth transition toward the depth cap |
| `exponential` | `D*(1-exp(-g*d/D))` | More gradual approach to the depth cap |

All start at zero depth at the bank with initial slope `g`. Maximum depth is a
ceiling, not a requirement: narrow channels remain shallower. The smooth profiles
approach their ceiling asymptotically. These describe the depth field relative
to the supplied water surface, not a bound on the final bed's total slope.
Tanh preserves the useful profile from the existing lake/coast code.

## Geometry and boundaries

- Supply valid, finite 2D Polygon or MultiPolygon geometry in the grid's projected
  metre-based CRS. No implicit reprojection or geometry repair occurs.
- Island holes remain land; their shores also influence depth.
- Grid dimensions, resolution and rotation are unrestricted. Sampling uses pixel
  centres; features narrower than a cell can disappear. No subpixel coverage or
  terrain smoothing is implied.
- Distances are measured directly against vector bank geometry in world metres,
  avoiding raster-distance stair steps and arbitrary smoothing constants.
- **Keep geometry beyond the grid edge.** The grid is a sampling window, not a
  shore. Clipping a polygon first creates false bank segments.
- For an open channel or artificially closed coastline, pass `banks=` containing
  only real shore segments as a LineString/MultiLineString. They must lie on the
  water boundary within 1 micrometre numerical tolerance. This prevents upstream,
  downstream or artificial closing edges from becoming raised sills. Include
  island shores yourself when overriding the default boundary.
- The water mask still requires a closed polygon. Constructing that polygon from
  an open coastline and deciding which side is water belong to the caller.

```python
from shapely.geometry import MultiLineString

channel = box(0, 0, 1000, 100)
real_banks = MultiLineString([
    [(0, 0), (1000, 0)],
    [(0, 100), (1000, 100)],
])
# generate_bed(channel, matching_grid, water_level_m=..., banks=real_banks)
# The two ends constrain coverage, but do not force depth to zero.
```

## Scope and extension points

This package generates deterministic **absolute target beds**, not successive
subtractions from a DEM. It does not yet apply or blend beds into terrain, detect
existing channels, estimate river surface elevations, reject steep reaches,
resolve overlapping water bodies, connect centrelines, place sources, or export
heightmaps. Lakes, rivers and coasts will be integrated separately later.

Geometry selection, reference-surface estimation, depth profiles and terrain
application are separate concerns. Future asymmetric channels or scour models
can extend depth generation without introducing data downloads or game-specific
logic here. No plugin framework is needed for the initial three profiles.

Calculations use float64; output requires about 17 bytes per cell for mask, depth
and elevation, plus temporary arrays. `chunk_rows` (default 128) bounds temporary
point geometry allocation, not output size. Exact vector distances favour simple,
accurate semantics; complex national-scale geometries may need future profiling
and spatial indexing. There are no performance guarantees for full map builds yet.

## Review plots

From the repository root, using the workspace's plotting dependencies:

```sh
uv run --all-packages python packages/water-carver/examples/plot_samples.py
```

Writes five annotated figures and their parameters to
`output/water-carver-samples/`: profile comparisons, an island lake, a winding
river, open versus closed channel ends, and the saved Hollyford polygon (when
available). The 3D views exaggerate vertical scale as labelled. Hollyford uses
an explicitly assumed flat reference surface, not terrain-derived water levels;
its mouth and lake ends remain untreated banks in this geometry-only example.

These plots also expose current limitations: wider reaches become deeper under
a nearest-bank rule, the nearest-bank junction can create a centre crease, and
a sloping reference surface does not ensure a monotonically descending bed.
They are bed-shape demonstrations, not hydraulic or in-game flow validation.


The selected experimental river settings are saved in
[`water-configs/river-preferred.json`](../../water-configs/river-preferred.json).
See the [river lab status](../../docs/river-lab.md#current-polygon-carver-experiment)
for the two-pass profile and Hollyford opening diagnostic. The profile and smoother are available in `water_carver.river`. The Jamestown
polygon trial uses these settings; the default generic `BedParameters` is unchanged.
