"""Hollyford parameter study. Run from repo root with uv run --all-packages."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from pyproj import Transformer
from shapely.affinity import translate
from shapely.geometry import shape
from shapely.ops import transform, unary_union
from scipy.ndimage import map_coordinates
from water_carver import BedParameters, generate_bed
from plot_samples import ROOT, grid_for, coordinates

OUT=ROOT/'output/water-carver-samples/hollyford'
GRADIENTS=[.05,.0633333333333333,.0866666666666667,.1]
CAPS=[5,10,15,20]
COLORS=['#2679af','#009688','#b77812','#b34769']
OUT.mkdir(parents=True,exist_ok=True)
source=ROOT/'river-configs/jamestown-hollyford-polygons.geojson'
data=json.loads(source.read_text())
project=Transformer.from_crs(4326,2193,always_xy=True).transform
geom=unary_union([transform(project,shape(f['geometry'])) for f in data['features']])
x0,y0,_,_=geom.bounds; geom=translate(geom,-x0,-y0)
_,_,x1,y1=geom.bounds
grid=grid_for((-30,-30,x1+30,y1+30),4)
beds={};stats=[]
for cap in CAPS:
    for g in GRADIENTS:
        bed=generate_bed(geom,grid,water_level_m=0,parameters=BedParameters(cap,g,'tanh'))
        beds[cap,g]=bed
        stats.append(dict(max_depth_cap_m=cap,bank_gradient=g,achieved_max_depth_m=float(np.nanmax(bed.depth_m)),mean_depth_m=float(np.nanmean(bed.depth_m))))


def plan(ax,bed,vmax):
    a=grid.transform;h,w=grid.shape
    return ax.imshow(bed.depth_m,extent=(a.c,a.c+4*w,a.f-4*h,a.f),origin='upper',interpolation='nearest',cmap='YlGnBu',vmin=0,vmax=vmax)

note='Actual LINZ polygon (CC BY 4.0); synthetic depths; flat reference water level 0 m; 4 m grid. All edges remain banks, including river ends.'
# Matrix uses one scale: no per-panel normalisation.
fig,axes=plt.subplots(4,4,figsize=(17,22),layout='constrained',sharex=True,sharey=True)
for r,cap in enumerate(CAPS):
    for c,g in enumerate(GRADIENTS):
        ax=axes[r,c];bed=beds[cap,g];im=plan(ax,bed,20)
        ax.set_title(f'{cap} m cap · {g*100:.2f}% bank\nDeepest cell: {np.nanmax(bed.depth_m):.2f} m',fontsize=10)
        if c==0:ax.set_ylabel('Local northing (m)')
        if r==3:ax.set_xlabel('Local easting (m)')
fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.35,label='Depth below assumed water surface (m) · same scale throughout')
fig.suptitle('Hollyford / tanh parameter study\nColumns: bank gradient     Rows: maximum depth cap',fontsize=20)
fig.text(.5,.005,note,ha='center',fontsize=9)
fig.set_constrained_layout_pads(h_pad=.4,w_pad=.1)
fig.savefig(OUT/'01-full-matrix.png',dpi=160);plt.close(fig)
# A larger 10m-cap comparison, useful at ordinary screen size.
fig,axes=plt.subplots(1,4,figsize=(15,11),layout='constrained',sharey=True)
for ax,g in zip(axes,GRADIENTS):
    im=plan(ax,beds[10,g],10)
    ax.set(title=f'{g*100:.2f}% bank gradient\nDeepest: {np.nanmax(beds[10,g].depth_m):.2f} m',xlabel='Local easting (m)')
axes[0].set_ylabel('Local northing (m)')
fig.colorbar(im,ax=list(axes),shrink=.5,label='Depth (m) · shared 0–10 m scale')
fig.suptitle('Hollyford / fixed 10 m cap, varying bank gradient',fontsize=19)
fig.text(.5,.005,note,ha='center',fontsize=8);fig.set_constrained_layout_pads(h_pad=.5,w_pad=.1)
fig.savefig(OUT/'02-ten-metre-cap.png',dpi=170);plt.close(fig)
# Exact analytic profiles. 250m bank distance approximately covers this polygon.
d=np.linspace(0,250,1001)
fig,axes=plt.subplots(2,4,figsize=(17,9),layout='constrained',sharex=True,sharey=True)
for c,cap in enumerate(CAPS):
    for g,color in zip(GRADIENTS,COLORS):
        axes[0,c].plot(d,-BedParameters(cap,g).depths(d),label=f'{g*100:.2f}%',color=color,lw=2)
    axes[0,c].set_title(f'Fixed cap: {cap} m');axes[0,c].legend(title='Bank gradient',fontsize=8)
for c,g in enumerate(GRADIENTS):
    for cap,color in zip(CAPS,COLORS):
        axes[1,c].plot(d,-BedParameters(cap,g).depths(d),label=f'{cap} m',color=color,lw=2)
    axes[1,c].set_title(f'Fixed bank gradient: {g*100:.2f}%');axes[1,c].legend(title='Depth cap',fontsize=8)
for ax in axes.flat:
    ax.set(xlabel='Distance inward from nearest bank (m)',ylabel='Bed relative to water surface (m)',ylim=(-21,.5));ax.grid(alpha=.2)
fig.suptitle('Tanh profiles / separate the effect of bank gradient from the depth cap',fontsize=18)
fig.savefig(OUT/'03-profile-comparisons.png',dpi=170);plt.close(fig)
# Fixed E–W transects in actual polygon; islands/gaps remain NaN.
sections=[('Broad northern reach',4850),('Island reach',4000),('Southern reach',1250)]
fig,axes=plt.subplots(3,3,figsize=(17,12),layout='constrained',width_ratios=[.8,1.3,1.3])
x=np.arange(-20,x1+20,2)
for r,(label,y) in enumerate(sections):
    im=plan(axes[r,0],beds[10,GRADIENTS[1]],10)
    axes[r,0].axhline(y,color='#ba4946',ls='--',lw=1.5)
    axes[r,0].set(title=label,xlabel='Easting (m)',ylabel='Northing (m)')
    mask=map_coordinates(beds[10,GRADIENTS[1]].mask.astype(float),[np.full_like(x,(grid.transform.f-y)/4-.5),(x-grid.transform.c)/4-.5],order=0,mode='constant',cval=0)>0
    valid=x[mask]
    for g,color in zip(GRADIENTS,COLORS):
        vals=map_coordinates(beds[10,g].elevation_m,[np.full_like(x,(grid.transform.f-y)/4-.5),(x-grid.transform.c)/4-.5],order=1,mode='constant',cval=np.nan)
        axes[r,1].plot(x,vals,label=f'{g*100:.2f}%',color=color,lw=2)
    for cap,color in zip(CAPS,COLORS):
        vals=map_coordinates(beds[cap,GRADIENTS[1]].elevation_m,[np.full_like(x,(grid.transform.f-y)/4-.5),(x-grid.transform.c)/4-.5],order=1,mode='constant',cval=np.nan)
        axes[r,2].plot(x,vals,label=f'{cap} m',color=color,lw=2)
    for ax,title in zip(axes[r,1:],['Fixed 10 m cap · vary gradient','Fixed 6.33% gradient · vary cap']):
        ax.axhline(0,color='#597585',ls='--',lw=1)
        ax.set(title=title,xlabel='Local easting (m)',ylabel='Bed elevation (m)',ylim=(-16,.7),xlim=(valid.min()-20,valid.max()+20));ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Actual Hollyford cross-sections / dashed lines locate the cuts\nEast–west sections, not width estimates; gaps preserve islands or land',fontsize=18)
fig.savefig(OUT/'04-cross-sections.png',dpi=160);plt.close(fig)
(OUT/'parameters.json').write_text(json.dumps(dict(profile='tanh',gradients=GRADIENTS,caps_m=CAPS,reference_level_m=0,resolution_m=4,assumptions=note,statistics=stats),indent=2)+'\n')
(OUT/'README.md').write_text('''# Hollyford tanh comparisons

All figures use the actual saved Hollyford polygon. Water level is assumed flat at 0 m; no terrain is modified. River-end banks have not been removed. These are synthetic shape experiments, not surveyed depths or water-flow predictions.

## Fixed 10 m cap

![Fixed 10 m cap](02-ten-metre-cap.png)

## All 16 combinations

The shared 0–20 m colour scale makes depths comparable across all panels. The enlarged 10 m-cap figure above instead uses a shared 0–10 m scale.

![Full matrix](01-full-matrix.png)

## Profile comparisons

![Profiles](03-profile-comparisons.png)

## Actual polygon cross-sections

![Cross-sections](04-cross-sections.png)

Depth caps are ceilings, not guaranteed depths. Increasing the cap has progressively less effect where the available distance to shore is small. Bank gradients are m/m percentages, not degrees or downstream slopes.
''')
print(json.dumps(stats,indent=2))
