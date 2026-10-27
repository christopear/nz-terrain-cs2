"""Endpoint-aligned cubic curve; no iterative optimisation or terrain changes."""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d
from shapely.geometry import LineString


def fit_curve(route, connection, anchor_depth=20., fit_radius=15., guide_spacing=40., guide_smoothing=10., output_spacing=1.):
    if min(anchor_depth,fit_radius,guide_spacing,output_spacing)<=0 or guide_smoothing<0:
        raise ValueError('Distances must be positive; smoothing may be zero')
    endpoints=[]
    for station,sign in zip(connection.stations,(-1,1)):
        # Group by proximity to the corresponding crossing, independent of input ID order.
        groups={s['end'] for s in connection.samples}
        label=min(groups,key=lambda label:min(abs(route.project(s['geometry'].interpolate(.5,normalized=True))-station) for s in connection.samples if s['end']==label))
        group=[s for s in connection.samples if s['end']==label]
        depth=min(anchor_depth,max(s['offset_m'] for s in group))
        selected=[s for s in group if abs(s['offset_m']-depth)<=fit_radius]
        if len(selected)<3: raise ValueError('Need at least three nearby transects for endpoint direction')
        x=np.array([sign*(s['offset_m']-depth) for s in selected])
        midpoints=np.array([s['geometry'].interpolate(.5,normalized=True).coords[0] for s in selected])
        slope,position=np.linalg.lstsq(np.column_stack([x,np.ones(len(x))]),midpoints,rcond=None)[0]
        if np.linalg.norm(slope)<1e-8: raise ValueError('Degenerate midpoint direction')
        direction=slope/np.linalg.norm(slope)
        endpoints.append((station+sign*depth,position,direction,midpoints))
    start,end=endpoints
    # Smooth the mapped route modestly, then retain guide points through the gap.
    stations=np.linspace(start[0],end[0],max(3,int((end[0]-start[0])/2)+1))
    xy=np.array([route.interpolate(s).coords[0] for s in stations])
    if guide_smoothing: xy=gaussian_filter1d(xy,guide_smoothing/(stations[1]-stations[0]),axis=0,mode='nearest')
    inner=np.arange(connection.stations[0]+guide_spacing,connection.stations[1]-guide_spacing/2,guide_spacing)
    inner=inner[(inner>start[0]) & (inner<end[0])]
    knots=np.r_[start[0],inner,end[0]]
    positions=np.vstack([start[1],np.column_stack([np.interp(inner,stations,xy[:,k]) for k in (0,1)]),end[1]])
    spline=CubicSpline(knots,positions,bc_type=((1,start[2]),(1,end[2])))
    dense=np.linspace(start[0],end[0],max(3,int(np.ceil((end[0]-start[0])/output_spacing))+1))
    curve=LineString(spline(dense))
    if not curve.is_simple: raise ValueError('Curve self-intersects; inspect guide spacing')
    corridor=curve.buffer(connection.width/2,cap_style='flat',join_style='round')
    return curve,corridor,endpoints,positions
