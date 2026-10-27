"""Metre-based river cross-sections with explicit unknown-width outcomes."""
import numpy as np
from scipy.ndimage import map_coordinates
from shapely.geometry import LineString, Point


def line_parts(geometry):
    if geometry.is_empty:
        return []
    if geometry.geom_type == 'LineString':
        return [geometry]
    return [line for part in getattr(geometry, 'geoms', []) for line in line_parts(part)]


def transect(line, distance, half_width, tangent_window=25.):
    """Construct a perpendicular section using a local, smoothed route tangent."""
    point = line.interpolate(distance)
    before = np.array(line.interpolate(max(0, distance-tangent_window)).coords[0][:2])
    after = np.array(line.interpolate(min(line.length, distance+tangent_window)).coords[0][:2])
    delta = after-before
    if np.linalg.norm(delta) < 1e-8:
        raise ValueError('Cannot construct a transect at a degenerate route tangent')
    normal = np.array([-delta[1], delta[0]]) / np.linalg.norm(delta)
    xy = np.array(point.coords[0][:2])
    return LineString([xy-half_width*normal, xy+half_width*normal])


def polygon_width(section, centre, polygon):
    """Measure only the connected polygon interval containing the route.

    Do not bridge islands, sum disconnected braids, or substitute zero for
    missing evidence. A section ending in water is censored, not a full width.
    """
    if not polygon.covers(centre):
        return None, 'centre_outside_polygon'
    pieces = [p for p in line_parts(section.intersection(polygon)) if p.distance(centre) < 1e-6]
    if len(pieces) != 1:
        return None, 'ambiguous_intersection'
    piece = pieces[0]
    if any(piece.distance(Point(xy)) < 1e-6 for xy in [section.coords[0], section.coords[-1]]):
        return None, 'transect_truncated'
    return float(piece.length), 'mapped_polygon'


def sample_surface(surface, xy, affine):
    """Bilinear sampling with NaN beyond the original DEM's pixel centres."""
    columns, rows = (~affine) * (xy[:, 0], xy[:, 1])
    return map_coordinates(surface, [rows-.5, columns-.5], order=1, mode='constant', cval=np.nan)
