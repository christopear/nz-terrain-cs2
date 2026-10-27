"""Boundary-distance bed profiles, independent of terrain and map pipelines."""
from dataclasses import dataclass
from numbers import Integral
from typing import Literal

from affine import Affine
import numpy as np
from numpy.typing import ArrayLike, NDArray
from pyproj import CRS
import shapely
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon


@dataclass(frozen=True)
class Grid:
    """Raster shape and pixel-corner affine in a projected, metre-based CRS.

    Polygon and bank coordinates must already be in this CRS. Rotated and
    sheared grids are supported; distances are measured in world coordinates.
    """

    shape: tuple[int, int]
    transform: Affine
    crs: CRS | str | int

    def __post_init__(self):
        if len(self.shape) != 2 or any(
            isinstance(n, bool) or not isinstance(n, Integral) or n <= 0
            for n in self.shape
        ):
            raise ValueError("Grid shape must contain two positive integers")
        if not isinstance(self.transform, Affine):
            raise ValueError("Grid transform must be an Affine")
        if not np.isfinite(tuple(self.transform)).all() or self.transform.determinant == 0:
            raise ValueError("Grid transform must be finite and invertible")
        crs = CRS.from_user_input(self.crs)
        if not crs.is_projected or len(crs.axis_info) < 2 or any(
            not np.isclose(axis.unit_conversion_factor, 1., rtol=0, atol=1e-12)
            for axis in crs.axis_info[:2]
        ):
            raise ValueError("Use a projected CRS with horizontal units in metres")
        object.__setattr__(self, "crs", crs)
        object.__setattr__(self, "shape", tuple(self.shape))


@dataclass(frozen=True)
class BedParameters:
    """Depth in metres; bank_gradient is vertical metres per horizontal metre.

    All profiles start at depth zero at a bank. Tanh and exponential approach
    max_depth_m asymptotically; linear reaches it and then stays flat.
    bank_gradient is the initial depth slope, not downstream river slope.
    """

    max_depth_m: float = 30.
    bank_gradient: float = .05
    profile: Literal["tanh", "linear", "exponential"] = "tanh"

    def __post_init__(self):
        if any(not np.isfinite(v) or v <= 0 for v in (self.max_depth_m, self.bank_gradient)):
            raise ValueError("Maximum depth and bank gradient must be finite and positive")
        if self.profile not in ("tanh", "linear", "exponential"):
            raise ValueError("Profile must be tanh, linear or exponential")

    def depths(self, distance_m: ArrayLike) -> NDArray[np.float64]:
        """Evaluate the profile for finite, nonnegative distances."""
        distance = np.asarray(distance_m, dtype=np.float64)
        if not np.isfinite(distance).all() or (distance < 0).any():
            raise ValueError("Bank distances must be finite and nonnegative")
        with np.errstate(over="ignore"):
            scaled = distance * self.bank_gradient / self.max_depth_m
        if self.profile == "linear":
            fraction = np.minimum(scaled, 1.)
        elif self.profile == "tanh":
            fraction = np.tanh(scaled)
        else:
            fraction = -np.expm1(-scaled)
        return self.max_depth_m * fraction


@dataclass(frozen=True)
class Bed:
    """Absolute elevations and positive depths; NaN outside the water mask.

    Float64 arrays retain precision when small depths are subtracted from high
    elevations. The caller owns export, blending and terrain application.
    """

    grid: Grid
    mask: NDArray[np.bool_]
    depth_m: NDArray[np.float64]
    elevation_m: NDArray[np.float64]


def generate_bed(
    water: Polygon | MultiPolygon,
    grid: Grid,
    *,
    water_level_m: ArrayLike,
    parameters: BedParameters | None = None,
    banks: LineString | MultiLineString | None = None,
    chunk_rows: int = 128,
) -> Bed:
    """Sample a polygon bed at raster pixel centres.

    Reference water levels may be a scalar or an array matching grid.shape.
    Values inside water must be finite and share the desired bed's vertical
    datum. No DEM inference, excavation, slope filtering or hydrology occurs.

    By default every polygon boundary (including island shores) is a bank.
    Supply only the true bank portions for open channels/artificial closures:
    omitted edges then constrain the mask but do not force depth back to zero.
    Banks must lie on the polygon boundary (1 micrometre numerical tolerance).
    Never clip water geometry to the raster before calling this function:
    doing so would introduce artificial banks at the map edge.
    """
    if not isinstance(water, (Polygon, MultiPolygon)) or water.is_empty or not water.is_valid:
        raise ValueError("Water must be a nonempty valid Polygon or MultiPolygon")
    if not np.isfinite(shapely.get_coordinates(water)).all() or shapely.has_z(water):
        raise ValueError("Water must have finite 2D coordinates")
    parameters = parameters if parameters is not None else BedParameters()
    boundary = water.boundary
    if banks is not None:
        if not isinstance(banks, (LineString, MultiLineString)) or banks.is_empty or not banks.is_valid:
            raise ValueError("Banks must be nonempty valid line geometry")
        if shapely.has_z(banks) or not np.isfinite(shapely.get_coordinates(banks)).all():
            raise ValueError("Banks must have finite 2D coordinates")
        if not boundary.buffer(1e-6).covers(banks):
            raise ValueError("Banks must lie on the water polygon boundary")
        boundary = banks
    if isinstance(chunk_rows, bool) or not isinstance(chunk_rows, Integral) or chunk_rows <= 0:
        raise ValueError("chunk_rows must be a positive integer")
    level = np.asarray(water_level_m, dtype=np.float64)
    if level.ndim != 0 and level.shape != grid.shape:
        raise ValueError("Water level must be scalar or match the grid shape")
    if level.ndim == 0 and not np.isfinite(level):
        raise ValueError("Scalar water level must be finite")
    mask = np.zeros(grid.shape, dtype=bool)
    depth = np.full(grid.shape, np.nan)
    elevation = np.full(grid.shape, np.nan)
    affine = grid.transform
    columns = np.arange(grid.shape[1], dtype=np.float64)[None, :] + .5
    # Vector distances avoid raster stair-step shores, padding heuristics and
    # resolution-dependent smoothing. Chunking bounds temporary geometry use.
    for start in range(0, grid.shape[0], chunk_rows):
        stop = min(start + chunk_rows, grid.shape[0])
        rows = np.arange(start, stop, dtype=np.float64)[:, None] + .5
        x = affine.a * columns + affine.b * rows + affine.c
        y = affine.d * columns + affine.e * rows + affine.f
        points = shapely.points(x, y)
        inside = shapely.covers(water, points)
        mask[start:stop] = inside
        if not inside.any():
            continue
        reference = level if level.ndim == 0 else level[start:stop][inside]
        if not np.isfinite(reference).all():
            raise ValueError("Water levels inside the water polygon must be finite")
        values = parameters.depths(shapely.distance(points[inside], boundary))
        depth[start:stop][inside] = values
        elevation[start:stop][inside] = reference - values
    return Bed(grid, mask, depth, elevation)
