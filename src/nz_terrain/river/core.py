"""Legacy numerical primitives, independent of presets and water-bed builders.

These retain the existing carving behaviour for compatibility. The Jamestown
lab measures inputs separately; it does not treat this carver as the new design.
"""
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.optimize import linprog
from scipy.sparse import lil_matrix


def solve_profile(reference, distance, depth, max_cut, slope, outlet, source=None):
    """Highest feasible cut-only bed with a soft minimum downhill slope.

    Reference is a smoothed surface proxy, not a measured water level. Fixed
    downstream bed and hard non-uphill constraints apply to the whole reach.
    Slack absorbs unavailable drop rather than excavating an unbounded trench.
    """
    n = len(reference)
    upper = reference-depth
    lower = reference-max_cut
    if source is not None:
        lower=np.maximum(lower,source-max_cut+.1)
        upper=np.minimum(upper,source-.25)
    lower[-1] = upper[-1] = outlet
    # Variables: bed elevations, followed by minimum-slope shortfalls (metres).
    a = lil_matrix((2*(n-1), 2*n-1))
    for i in range(n-1):
        a[i,i], a[i,i+1] = -1, 1
        a[n-1+i,i], a[n-1+i,i+1], a[n-1+i,n+i] = -1, 1, -1
    drop = np.diff(distance)*slope
    # Maximise bed elevation; a modest penalty favours the target slope without
    # sacrificing metres of extra excavation to meet it over a long flat reach.
    cost = np.r_[-np.ones(n), np.full(n-1, 20.)]
    result = linprog(cost, A_ub=a.tocsr(), b_ub=np.r_[np.zeros(n-1), -drop],
                     bounds=list(zip(lower,upper))+[(0,None)]*(n-1),method='highs')
    if not result.success:
        raise ValueError('River profile cannot meet depth/cut/outlet constraints: '+result.message)
    bed = result.x[:n]
    assert np.max(np.diff(bed)) < 1e-6
    return bed


def sample(z,xy,affine):
    c,r=(~affine)*(xy[:,0],xy[:,1])
    return map_coordinates(z,[r-.5,c-.5],order=1,mode='nearest')


def carve(z,profiles,config,affine):
    """Evaluate continuous segment distance/profiles in bounded pixel windows."""
    result=z.copy();mask=np.zeros(z.shape,bool)
    for name,p in profiles.items():
        dims=config['rivers'][name]
        for i in range(len(p['xy'])-1):
            bank=p.get('bank_width',np.full(len(p['xy']),dims['bank_width_m']))[i:i+2]
            bottom=p.get('bottom_width',np.full(len(p['xy']),dims['bottom_width_m']))[i:i+2]
            half=float(bank.max())/2
            a,b=p['xy'][i:i+2];delta=b-a;length2=delta@delta
            if length2<1e-8:continue
            cc,rr=(~affine)*(np.array([a[0],b[0]]),np.array([a[1],b[1]]))
            pad=half/3.5+2
            c0=max(0,int(np.floor(cc.min()-pad)));c1=min(z.shape[1],int(np.ceil(cc.max()+pad)))
            r0=max(0,int(np.floor(rr.min()-pad)));r1=min(z.shape[0],int(np.ceil(rr.max()+pad)))
            yy,xx=np.mgrid[r0:r1,c0:c1];x,y=affine*(xx+.5,yy+.5)
            t=np.clip(((x-a[0])*delta[0]+(y-a[1])*delta[1])/length2,0,1)
            half=(bank[0]*(1-t)+bank[1]*t)/2
            flat=(bottom[0]*(1-t)+bottom[1]*t)/2
            dist=np.hypot(x-a[0]-t*delta[0],y-a[1]-t*delta[1])
            inside=dist<half
            bed=p['bed'][i]*(1-t)+p['bed'][i+1]*t
            # Flat bottom and cosine bank blend meet unchanged terrain smoothly.
            blend=(1+np.cos(np.pi*np.clip((dist-flat)/(half-flat),0,1)))/2
            original=z[r0:r1,c0:c1]
            target=original+blend*(np.minimum(original,bed)-original)
            # Cap excavation at each raster pixel, including steep banks. Any
            # resulting centreline obstruction is measured in the report.
            target=np.maximum(target,original-dims.get('max_cut_m',config['max_cut_m']))
            view=result[r0:r1,c0:c1]
            view[inside]=np.minimum(view[inside],target[inside]);mask[r0:r1,c0:c1]|=inside
    return result,mask
