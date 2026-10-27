"""Bank-driven river depth and local masked crease diffusion, in metre units."""
import numpy as np
from scipy.ndimage import distance_transform_edt

def profile(distance,transition,cap=10,shore_gradient=.05,deep_gradient=.5):
    """S-like depth: shallow initial slope, steepening bank, then tanh cap.

    Bank slope before cap = shore_gradient + (deep_gradient-shore_gradient)
    * smoothstep(distance/transition). Its integral is used as the tanh input.
    Therefore transition is an exact horizontal distance, not a sigmoid scale.
    """
    u=np.clip(distance/transition,0,1)
    integral=transition*(u**3-.5*u**4)+np.maximum(distance-transition,0)
    descent=shore_gradient*distance+(deep_gradient-shore_gradient)*integral
    return cap*np.tanh(descent/cap)


def smooth_creases(depth,distance,spacing,sigma,shore_lock=.5,shore_blend=1.5):
    """Experimental local, masked diffusion. Does not read channel widths.

    Find distance-field ridges via negative discrete Laplacian; diffuse depth
    only across adjacent water cells, locking near-shore samples. Fade diffusion
    strength to zero at the edge of a finite 2*sigma band around ridges. Returns depth and blend strength.
    Smoothing fills sharp deepest points; it does not promise to preserve depth
    or volume. Frozen land never transmits diffusion between channels.
    """
    water=np.isfinite(depth)
    d=np.where(water,distance,0.)
    lap=np.zeros_like(d)
    for axis in range(d.ndim):
        lo=[slice(None)]*d.ndim;hi=lo.copy();lo[axis]=slice(None,-1);hi[axis]=slice(1,None)
        lo=tuple(lo);hi=tuple(hi);pair=water[lo]&water[hi]
        delta=np.where(pair,d[hi]-d[lo],0.)
        lap[lo]+=delta;lap[hi]-=delta
    ridge=water & (lap < -.1*spacing) & (d>shore_lock)
    if not ridge.any():return depth.copy(),np.zeros_like(d)
    proximity=distance_transform_edt(~ridge)*spacing
    u=np.clip(1-proximity/(2*sigma),0,1)
    ridge_weight=u*u*(3-2*u)
    u=np.clip((d-shore_lock)/shore_blend,0,1)
    weight=ridge_weight*u*u*(3-2*u)*water
    current=np.where(water,depth,0.)
    active=water & (d>shore_lock)
    # alpha <= 1/(2*ndim) gives a convex, stable explicit diffusion update.
    total_time=sigma*sigma/(2*spacing*spacing)
    steps=max(1,int(np.ceil(total_time/.2)))
    alpha=total_time/steps
    for _ in range(steps):
        flux=np.zeros_like(current)
        for axis in range(current.ndim):
            lo=[slice(None)]*current.ndim;hi=lo.copy();lo[axis]=slice(None,-1);hi[axis]=slice(1,None)
            lo=tuple(lo);hi=tuple(hi);pair=water[lo]&water[hi]
            delta=np.where(pair,current[hi]-current[lo],0.)
            flux[lo]+=delta;flux[hi]-=delta
        current[active]+=alpha*weight[active]*flux[active]
    result=np.where(water,current,np.nan)
    return result,weight
