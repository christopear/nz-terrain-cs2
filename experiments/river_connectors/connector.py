"""Standalone geometry experiment. Coordinates and parameters must be in metres."""
from dataclasses import dataclass
import numpy as np
from shapely.ops import substring
from nz_terrain.river.geometry import line_parts, transect, polygon_width


@dataclass
class Connection:
    gap: object
    corridor: object
    samples: list
    end_widths: tuple
    width: float
    stations: tuple


def connect(route, first, second, sample_run=75., sample_step=5., tangent_window=10.):
    """Bridge the shortest along-route gap between two polygons.

    Average perpendicular sections inside each adjacent polygon reach, then give
    each polygon equal weight in the constant connector width. No DEM or carving.
    A short overlap into each polygon makes the proposed corridor joinable.
    """
    if route.geom_type != 'LineString' or not route.is_simple:
        raise ValueError('Requires one simple, continuous centreline')
    if min(sample_run, sample_step, tangent_window) <= 0:
        raise ValueError('Sampling parameters must be positive')
    intervals=[]
    for label,poly in enumerate((first,second)):
        for segment in line_parts(route.intersection(poly)):
            a,b=sorted(route.project(segment.interpolate(s)) for s in (0,segment.length))
            if b-a>1e-6: intervals.append((a,b,label))
    intervals.sort()
    candidates=[(left,right) for left,right in zip(intervals,intervals[1:])
                if left[2]!=right[2] and right[0]-left[1]>1e-6]
    if not candidates: raise ValueError('No uncovered route gap between these polygons')
    left,right=min(candidates,key=lambda pair:pair[1][0]-pair[0][1])
    samples=[];means=[]
    for interval,station,direction in ((left,left[1],-1),(right,right[0],1)):
        poly=(first,second)[interval[2]]
        widths=[]
        available=min(sample_run,interval[1]-interval[0])
        for offset in np.arange(sample_step/2,available,sample_step):
            s=station+direction*offset
            section=transect(route,s,1000.,tangent_window)
            width,status=polygon_width(section,route.interpolate(s),poly)
            if width is None: continue
            piece=next(p for p in line_parts(section.intersection(poly)) if p.distance(route.interpolate(s))<1e-6)
            widths.append(width)
            samples.append(dict(end=interval[2],offset_m=float(offset),width_m=width,geometry=piece))
        if not widths: raise ValueError('No usable width measurements inside a polygon')
        means.append(float(np.mean(widths)))
    width=float(np.mean(means)); start,end=left[1],right[0]
    # Overlap is limited to the measured adjacent polygon reaches.
    extended=substring(route,max(left[0],start-width),min(right[1],end+width))
    corridor=extended.buffer(width/2,cap_style='flat',join_style='round')
    return Connection(substring(route,start,end),corridor,samples,tuple(means),width,(start,end))
