"""Longer bank runs and a width-aware hybrid prototype; no production changes.
Run with: uv run --all-packages python packages/water-carver/examples/plot_long_banks.py
"""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer
from scipy.ndimage import map_coordinates
from shapely.affinity import translate
from shapely.geometry import shape
from shapely.ops import transform, unary_union
from water_carver import Bed, BedParameters, generate_bed
from plot_samples import ROOT, grid_for

OUT=ROOT/'output/water-carver-samples/hollyford-long-banks'
OUT.mkdir(parents=True,exist_ok=True)
project=Transformer.from_crs(4326,2193,always_xy=True).transform
source=ROOT/'river-configs/jamestown-hollyford-polygons.geojson'
water=unary_union([transform(project,shape(f['geometry'])) for f in json.loads(source.read_text())['features']])
x0,y0,_,_=water.bounds;water=translate(water,-x0,-y0)
_,_,x1,y1=water.bounds;grid=grid_for((-30,-30,x1+30,y1+30),2)
base=generate_bed(water,grid,water_level_m=0,parameters=BedParameters(np.hypot(x1,y1)+100,1,'linear'))
runs=[40,60,80,100];colors=['#7a52a4','#008778','#be7c20','#bf4d72']

def shelf(d,run):
    u=np.clip(d/run,0,1)
    return 10*u*u*(3-2*u)


def hybrid(d,width,flat_width,depth=10,nominal_gradient=.3):
    """1D only: width sets bank-run length, a Hermite tail flattens tanh.

    Tanh is normalised to the target depth at the floor edge; the last portion
    transitions to zero derivative using a cubic Hermite curve. This prioritises
    10m depth, so actual initial slope exceeds nominal_gradient in narrow cuts.
    Floor width is capped at one third of the channel width for tiny channels.
    """
    floor=min(flat_width,width/3)
    length=(width-floor)/2
    scale=nominal_gradient/depth
    norm=np.tanh(scale*length)
    blend=min(.4*length,depth/nominal_gradient)
    start=length-blend
    y0=depth*np.tanh(scale*start)/norm
    m0=nominal_gradient/(np.cosh(scale*start)**2*norm)
    u=np.clip((d-start)/blend,0,1)
    tail=(2*u**3-3*u**2+1)*y0+(u**3-2*u**2+u)*blend*m0+(-2*u**3+3*u**2)*depth
    values=np.where(d<start,depth*np.tanh(scale*d)/norm,tail)
    return np.where(d>=length,depth,values)

# Validate the experimental 1D construction, including small and wide channels.
for width in [10,20,40,60,100,200,400,1000]:
    for flat in [10,20]:
        d=np.linspace(0,width/2,2001);z=hybrid(d,width,flat)
        assert abs(z[0])<1e-12 and abs(z[-1]-10)<1e-12
        assert np.all(np.diff(z)>=-1e-10) and np.all(z<=10+1e-10)

beds={};stats=[]
for run in runs:
    z=np.full(grid.shape,np.nan);z[base.mask]=shelf(base.depth_m[base.mask],run)
    beds[run]=Bed(grid,base.mask,z,-z)
    stats.append(dict(bank_run_m=run,max_depth_m=float(np.nanmax(z)),mean_depth_m=float(np.nanmean(z)),max_bank_slope_percent=1500/run))

def plan(ax,bed):
    a=grid.transform;h,w=grid.shape
    return ax.imshow(bed.depth_m,origin='upper',extent=(a.c,a.c+2*w,a.f-2*h,a.f),vmin=0,vmax=10,cmap='YlGnBu',interpolation='nearest')

fig,axes=plt.subplots(1,4,figsize=(15,11),layout='constrained',sharey=True)
for ax,run in zip(axes,runs):
    im=plan(ax,beds[run]);ax.set(title=f'{run} m bank run\nPeak bank slope: {1500/run:g}%',xlabel='Easting (m)')
axes[0].set_ylabel('Northing (m)');fig.colorbar(im,ax=list(axes),shrink=.5,label='Synthetic depth (m) · shared scale')
fig.suptitle('Hollyford / rounded shelves, longer bank transitions\n10 m cap throughout · 40 m is the previous purple reference',fontsize=18)
fig.text(.5,.005,'LINZ polygon: CC BY 4.0. 2 m grid, assumed flat water level 0 m, all boundary edges remain banks. No terrain changes.',ha='center',fontsize=9)
fig.set_constrained_layout_pads(h_pad=.4,w_pad=.1)
fig.savefig(OUT/'01-longer-bank-polygons.png',dpi=170);plt.close(fig)

fig,axes=plt.subplots(3,2,figsize=(14,12),layout='constrained',width_ratios=[.6,1.4])
x=np.arange(-20,x1+20,1.)
for row,(title,y) in enumerate([('Broad northern reach',4850),('Island reach',4000),('Southern reach',1250)]):
    plan(axes[row,0],beds[60]);axes[row,0].axhline(y,color='#bd4141',ls='--');axes[row,0].set(title=title,xlabel='Easting (m)',ylabel='Northing (m)')
    coords=[np.full_like(x,(grid.transform.f-y)/2-.5),(x-grid.transform.c)/2-.5]
    mask=map_coordinates(base.mask.astype(float),coords,order=0,mode='constant',cval=0)>0
    valid=x[mask]
    for run,color in zip(runs,colors):
        z=map_coordinates(beds[run].elevation_m,coords,order=1,mode='constant',cval=np.nan)
        axes[row,1].plot(x,z,color=color,lw=2,label=f'{run} m bank run')
    ax=axes[row,1];ax.set(xlabel='Local easting (m)',ylabel='Bed elevation (m)',ylim=(-11,.5),xlim=(valid.min()-15,valid.max()+15));ax.legend(fontsize=9);ax.grid(alpha=.2)
fig.suptitle('Actual Hollyford cross-sections / 10 m cap\nLonger banks reduce depth where the channel cannot fit both transitions',fontsize=18)
fig.savefig(OUT/'02-actual-cross-sections.png',dpi=170);plt.close(fig)

fig,axes=plt.subplots(2,3,figsize=(17,10),layout='constrained',sharey=True)
for ax,width in zip(axes.flat,[40,60,100,160,240,400]):
    x=np.linspace(0,width,2001);d=np.minimum(x,width-x)
    ax.plot(x,-BedParameters(10,.3).depths(d),label='Original tanh · 10 m / 30%',color='#2b78b7',lw=2)
    ax.plot(x,-shelf(d,40),label='Previous rounded shelf · 40 m run',color=colors[0],lw=2)
    ax.plot(x,-hybrid(d,width,10),label='Hybrid · 10 m flat floor*',color='#008778',lw=2)
    ax.plot(x,-hybrid(d,width,20),label='Hybrid · 20 m flat floor*',color='#bd7023',lw=2,ls='--')
    ax.set(title=f'{width} m-wide ideal channel',xlabel='Across channel (m)',ylabel='Bed elevation (m)',ylim=(-11,.5));ax.grid(alpha=.2)
handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=2)
fig.suptitle('Width-aware tanh + rounded floor / CROSS-SECTION PROTOTYPE ONLY\n10 m target depth; final tanh segment eases into a flat bed',fontsize=18)
fig.text(.5,.08,'*Floor width limited to one third of channel width when necessary. Narrow channels require steeper-than-30% banks.\nThese ideal profiles have NOT been applied to the Hollyford polygon; local width estimation is still to be implemented.',ha='center',fontsize=10)
# Reserve a clear footer strip for caveats and legend.
fig.get_layout_engine().set(rect=(0,.12,1,.88))
fig.savefig(OUT/'03-width-aware-hybrid.png',dpi=170);plt.close(fig)
(OUT/'parameters.json').write_text(json.dumps(dict(fixed_run_variants=stats,depth_cap_m=10,hybrid_status='1D ideal-channel prototype only; not applied to polygon',hybrid_flat_widths_m=[10,20],hybrid_tiny_channel_rule='flat width <= channel width / 3',hybrid_nominal_gradient=.3,reference_water_level_m=0),indent=2)+'\n')
(OUT/'README.md').write_text('''# Longer banks and a width-aware hybrid

## Actual Hollyford: 40 / 60 / 80 / 100 m bank runs

The previous purple 40 m rounded shelf is the reference. Every variant has a 10 m cap. Longer bank runs make narrower reaches shallower and can reintroduce a central V when the two transitions meet before reaching the floor.

![Actual sections](02-actual-cross-sections.png)

![Polygon comparisons](01-longer-bank-polygons.png)

## Width-aware tanh + flat floor: 1D prototype only

![Hybrid profiles](03-width-aware-hybrid.png)

The hybrid allocates 10 or 20 m to a central flat floor, uses the remaining width for banks, and rounds the final part of tanh into that floor with a cubic Hermite transition. It retains a 10 m target depth. Tanh is normalised, so 30% is a nominal starting shape, NOT a guaranteed bank-slope limit. Narrow channels must get steeper to retain depth. Floor width is limited to one third of channel width for tiny sections.

This proves the cross-section shape only: no local-width map or 2D hybrid bed has been built yet. Applying it to actual polygons needs a consistent local-width estimate and smooth transitions at bends, junctions and islands. Fixed-distance shelves above do not need that machinery.

Actual polygon examples: LINZ CC BY 4.0, 2 m grid, flat assumed 0 m water surface, all polygon edges still act as banks. No DEM or water-flow simulation, and no production defaults changed.

Reproduce from repo root:

```sh
uv run --all-packages python packages/water-carver/examples/plot_long_banks.py
```
''')
print(json.dumps(stats,indent=2))
