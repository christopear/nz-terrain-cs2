"""Addition-only seam smoothing for synthetic connector polygons (metres)."""
from dataclasses import dataclass
from shapely.ops import unary_union


@dataclass
class JoinResult:
    original: object
    before: object
    additions: object
    synthetic_patch: object
    final: object
    seam_zone: object


def smooth_joins(polygons, connector, radius=5., seam_reach=25.):
    """Round concave seam notches, preserving all existing water geometry.

    Morphological closing (outward then inward round buffers) proposes fills.
    Only additions within seam_reach of source/connector boundary crossings
    are accepted. Source and connector geometry are never subtracted.
    This deliberately cannot round off protruding corners by removing area.
    """
    if radius<=0 or seam_reach<=radius:
        raise ValueError('Require 0 < radius < seam_reach')
    original=unary_union(polygons)
    before=original.union(connector)
    seams=original.boundary.intersection(connector.boundary)
    zone=seams.buffer(seam_reach,quad_segs=32)
    closed=before.buffer(radius,quad_segs=32).buffer(-radius,quad_segs=32)
    additions=closed.difference(before).intersection(zone)
    final=before.union(additions)
    patch=final.difference(original)
    return JoinResult(original,before,additions,patch,final,zone)
