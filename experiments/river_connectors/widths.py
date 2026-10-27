"""Variable-width corridor along an already fitted curve, in metres."""
import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union


def tapered_corridor(curve, end_widths, transition=200., spacing=1., mode="linear"):
    """Default to linear width interpolation along the entire curve.

    Optional end-transitions mode uses smoothstep ramps to the mean.

    Distances are curve arc length. On short reaches, each transition is limited
    to half the length so the two ramps meet at the mean without overlapping.
    Endpoint widths belong to curve order, not polygon argument order.
    """
    if min(*end_widths,transition,spacing,curve.length)<=0:
        raise ValueError('Widths, distances and curve length must be positive')
    if mode not in ("end-transitions", "linear"):
        raise ValueError("Unknown width interpolation mode")
    run=min(transition,curve.length/2)
    distances=np.linspace(0,curve.length,max(3,int(np.ceil(curve.length/spacing))+1))
    mid=float(np.mean(end_widths))
    def ramp(x):
        x=np.clip(x,0,1)
        return x*x*(3-2*x)
    widths=mid+(end_widths[0]-mid)*(1-ramp(distances/run))+(end_widths[1]-mid)*(1-ramp((curve.length-distances)/run))
    if mode=="linear":
        widths=np.interp(distances,[0,curve.length],end_widths)
        run=curve.length
    points=np.array([curve.interpolate(s).coords[0] for s in distances])
    before=np.array([curve.interpolate(max(0,s-.5)).coords[0] for s in distances])
    after=np.array([curve.interpolate(min(curve.length,s+.5)).coords[0] for s in distances])
    tangent=after-before
    norms=np.linalg.norm(tangent,axis=1)
    if np.any(norms<1e-9):raise ValueError('Degenerate curve tangent')
    normals=np.column_stack([-tangent[:,1],tangent[:,0]])/norms[:,None]
    left=points+normals*widths[:,None]/2
    right=points-normals*widths[:,None]/2
    # Union the swept strip: offset banks may fold over at tight inside bends,
    # just as a conventional buffer does. Triangles avoid invalid bow-tie quads.
    corridor=unary_union([triangle for i in range(len(left)-1) for triangle in (
        Polygon([left[i],right[i],left[i+1]]),
        Polygon([right[i],right[i+1],left[i+1]]))])
    if not corridor.is_valid or corridor.geom_type!='Polygon':raise ValueError('Swept corridor is not one valid polygon')
    return corridor,distances,widths,run
