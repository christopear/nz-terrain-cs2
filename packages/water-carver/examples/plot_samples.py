"""Render actual water-carver outputs; run from the repository root.

uv sync --all-packages --locked
uv run --all-packages python packages/water-carver/examples/plot_samples.py

Matplotlib is supplied by the nz-terrain workspace package. No source terrain
is read or modified. Outputs default to output/water-carver-samples/.
"""
import argparse
import json
from pathlib import Path

from affine import Affine
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer
from shapely.affinity import scale, translate
from scipy.ndimage import map_coordinates
from shapely.geometry import MultiLineString, Point, Polygon, box, shape
from shapely.ops import transform, unary_union

from water_carver import BedParameters, Grid, generate_bed

ROOT = Path(__file__).resolve().parents[3]
COLORS = {'linear':'#bd653b', 'tanh':'#007d88', 'exponential':'#7754ab'}
plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False,
                     'figure.facecolor':'#fcfcf8', 'axes.facecolor':'#f2f1eb',
                     'axes.titleweight':'bold', 'savefig.facecolor':'#fcfcf8'})


def grid_for(bounds, spacing=2):
    x0,y0,x1,y1=bounds
    rows=int(np.ceil((y1-y0)/spacing));cols=int(np.ceil((x1-x0)/spacing))
    return Grid((rows,cols),Affine(spacing,0,x0,0,-spacing,y1),'EPSG:2193')


def coordinates(grid):
    rows,cols=np.indices(grid.shape)
    a=grid.transform
    return a.a*(cols+.5)+a.b*(rows+.5)+a.c, a.d*(cols+.5)+a.e*(rows+.5)+a.f


def map_depth(ax, bed, maximum, title):
    a=bed.grid.transform;h,w=bed.grid.shape
    extent=(a.c,a.c+a.a*w,a.f+a.e*h,a.f)
    im=ax.imshow(bed.depth_m,extent=extent,origin='upper',cmap='YlGnBu',vmin=0,vmax=maximum,interpolation='nearest')
    ax.set(title=title,xlabel='Local easting (m)',ylabel='Local northing (m)',aspect='equal')
    return im


def save(fig,out,name):
    fig.savefig(out/name,dpi=170)
    plt.close(fig)


def render(out):
    out.mkdir(parents=True,exist_ok=True)
    metadata={'purpose':'Synthetic bed review; no DEM changes or water simulation', 'examples':{}}
    # 1. Profile behaviour and sensitivity, calculated by the production profile API.
    fig,axes=plt.subplots(1,3,figsize=(16,5),layout='constrained')
    d=np.linspace(0,800,801)
    for profile,color in COLORS.items():
        axes[0].plot(d,-BedParameters(30,.1,profile).depths(d),label=profile,color=color,lw=2.5)
    for g in [.025,.05,.1,.2]:
        axes[1].plot(d,-BedParameters(30,g).depths(d),label=f'{g*100:g}% bank gradient',lw=2)
    for cap in [5,15,30]:
        axes[2].plot(d,-BedParameters(cap,.1).depths(d),label=f'{cap} m depth cap',lw=2)
    for ax,title in zip(axes,['Profile shape · cap 30 m · gradient 10%', 'Bank gradient · tanh · cap 30 m','Depth cap · tanh · gradient 10%']):
        ax.set(title=title,xlabel='Distance inward from nearest bank (m)',ylabel='Bed relative to water surface (m)',ylim=(-32,1))
        ax.axhline(0,color='#597585',ls='--',lw=1);ax.grid(alpha=.2);ax.legend(loc='lower left',fontsize=9)
    fig.suptitle('01 / The same boundary, different bed profiles',fontsize=18)
    save(fig,out,'01-profiles.png')

    # 2. Lake with an island; holes remain unmodified/non-water.
    t=np.linspace(0,2*np.pi,361)
    outer=Polygon(np.column_stack([500+(430+28*np.sin(3*t))*np.cos(t),350+(275+18*np.cos(5*t))*np.sin(t)]))
    island=translate(scale(Point(0,0).buffer(1,quad_segs=48),65,43),380,365)
    lake=outer.difference(island);grid=grid_for((0,0,1000,700))
    lakebed=generate_bed(lake,grid,water_level_m=0,parameters=BedParameters(30,.12))
    fig=plt.figure(figsize=(15,9),layout='constrained');gs=fig.add_gridspec(2,2)
    ax=fig.add_subplot(gs[:,0]);im=map_depth(ax,lakebed,30,'Plan view · island remains land')
    ax.axhline(365,color='#b1483f',ls='--',lw=1.2);fig.colorbar(im,ax=ax,shrink=.6,label='Synthetic depth below water (m)')
    ax3=fig.add_subplot(gs[0,1],projection='3d');x,y=coordinates(grid)
    ax3.plot_surface(x[::4,::4],y[::4,::4],lakebed.elevation_m[::4,::4],cmap='terrain',linewidth=0,antialiased=True,rcount=200,ccount=200)
    ax3.set(xlabel='Easting (m)',ylabel='Northing (m)',zlabel='Bed (m)',title='3D bed · vertical scale exaggerated ~10×')
    ax3.set_zlim(-30,0);ax3.set_box_aspect((1000,700,300));ax3.view_init(32,-65)
    ax=fig.add_subplot(gs[1,1]);r=int((grid.transform.f-365)/2)
    ax.plot(x[r],lakebed.elevation_m[r],color=COLORS['tanh'],lw=2)
    ax.axhline(0,color='#597585',ls='--',label='Reference water level');ax.legend(fontsize=9)
    ax.set(title='Cross-section at dashed line · gap = island',xlabel='Local easting (m)',ylabel='Bed elevation (m)');ax.grid(alpha=.2)
    fig.suptitle('02 / Synthetic lake with an island\nTanh · depth cap 30 m · bank gradient 12% · water level 0 m',fontsize=16)
    save(fig,out,'02-lake-island.png')
    metadata['examples']['lake']={'max_depth_m':30,'bank_gradient':.12,'observed_max_depth_m':float(np.nanmax(lakebed.depth_m)),'resolution_m':2}

    # 3. Winding river polygon and explicit real banks; varying reference surface.
    xx=np.linspace(0,1200,241);centre=250+80*np.sin(xx/230);half=55+30*(1+np.sin(xx/170))/2
    lower=np.column_stack([xx,centre-half]);upper=np.column_stack([xx,centre+half])
    river=Polygon(np.concatenate([lower,upper[::-1]]));banks=MultiLineString([lower,upper])
    grid=grid_for((-20,0,1220,500));x,y=coordinates(grid);level=2-x*.001
    riverbed=generate_bed(river,grid,water_level_m=level,parameters=BedParameters(8,.08),banks=banks)
    fig=plt.figure(figsize=(15,9),layout='constrained');gs=fig.add_gridspec(2,2)
    ax=fig.add_subplot(gs[0,:]);im=map_depth(ax,riverbed,8,'Depth follows the actual polygon footprint · end caps excluded from banks')
    fig.colorbar(im,ax=ax,shrink=.6,label='Depth below local water surface (m)')
    ax=fig.add_subplot(gs[1,0]);dist=np.linspace(0,1200,601)
    route_y=250+80*np.sin(dist/230)
    # Display path only, not an input centreline: sample generated bed along it.
    cols=(dist-grid.transform.c)/2-.5
    rows=(grid.transform.f-route_y)/2-.5
    sampled_bed=map_coordinates(riverbed.elevation_m,[rows,cols],order=1,mode='constant',cval=np.nan)
    ax.plot(dist,2-dist*.001,color='#597585',ls='--',label='Supplied reference surface (0.1% fall)')
    ax.plot(dist,sampled_bed,color=COLORS['tanh'],label='Generated bed along display path',lw=2)
    ax.set(title='Along-channel profile · no flow solver',xlabel='Local easting (m)',ylabel='Elevation (m)');ax.legend(fontsize=8);ax.grid(alpha=.2)
    ax=fig.add_subplot(gs[1,1],projection='3d');stride=3
    ax.plot_surface(x[::stride,::stride],y[::stride,::stride],riverbed.elevation_m[::stride,::stride],cmap='terrain',linewidth=0,rcount=200,ccount=220)
    ax.set(xlabel='Easting (m)',ylabel='Northing (m)',zlabel='Bed (m)',title='3D bed · vertical scale exaggerated ~25×');ax.set_zlim(-6,2);ax.set_box_aspect((1240,500,200));ax.view_init(36,-65)
    fig.suptitle('03 / Synthetic river · varying width and sloping reference surface\nTanh · depth cap 8 m · bank gradient 8% · geometry only',fontsize=16)
    save(fig,out,'03-river.png')
    metadata['examples']['river']={'max_depth_m':8,'bank_gradient':.08,'reference_fall_m_per_m':.001,'observed_max_depth_m':float(np.nanmax(riverbed.depth_m)),'resolution_m':2}

    # 4. End caps are important: same rectangle, only bank selection differs.
    water=box(0,0,600,160);grid=grid_for((0,0,600,160));p=BedParameters(8,.1)
    closed=generate_bed(water,grid,water_level_m=0,parameters=p)
    opened=generate_bed(water,grid,water_level_m=0,parameters=p,banks=MultiLineString([[(0,0),(600,0)],[(0,160),(600,160)]]))
    fig,axes=plt.subplots(3,1,figsize=(12,9),layout='constrained')
    for ax,bed,title in zip(axes[:2],[closed,opened],['All edges treated as banks → bed rises across both ends','Only side banks → depth continues through both ends']):
        im=map_depth(ax,bed,8,title);fig.colorbar(im,ax=ax,shrink=.7,label='Depth (m)')
    x,_=coordinates(grid)
    for bed,label,color in [(closed,'Closed ends','#bd653b'),(opened,'Open ends','#007d88')]:
        axes[2].plot(x[40],bed.elevation_m[40],label=label,color=color,lw=2.5)
    axes[2].set(xlabel='Along channel (m)',ylabel='Centre bed elevation (m)',title='Longitudinal section · same polygon and depth parameters');axes[2].legend();axes[2].grid(alpha=.2)
    fig.suptitle('04 / Why open river ends must not be counted as banks\nTanh · depth cap 8 m · bank gradient 10% · water level 0 m',fontsize=16)
    save(fig,out,'04-open-ends.png')

    # 5. Actual saved Hollyford geometry, no terrain or inferred water levels.
    source=ROOT/'river-configs/jamestown-hollyford-polygons.geojson'
    if source.exists():
        data=json.loads(source.read_text());project=Transformer.from_crs(4326,2193,always_xy=True).transform
        geom=unary_union([transform(project,shape(f['geometry'])) for f in data['features']])
        x0,y0,_,_=geom.bounds;geom=translate(geom,-x0,-y0)
        x0,y0,x1,y1=geom.bounds;grid=grid_for((x0-30,y0-30,x1+30,y1+30),4)
        fig,axes=plt.subplots(1,3,figsize=(13,12),layout='constrained')
        stats=[]
        for ax,g in zip(axes,[.025,.05,.1]):
            bed=generate_bed(geom,grid,water_level_m=0,parameters=BedParameters(10,g))
            im=map_depth(ax,bed,10,f'Bank gradient {g*100:g}%\nActual max depth {np.nanmax(bed.depth_m):.1f} m')
            stats.append({'bank_gradient':g,'observed_max_depth_m':float(np.nanmax(bed.depth_m))})
        fig.colorbar(im,ax=list(axes),shrink=.5,label='Synthetic depth below assumed water surface (m)')
        fig.suptitle('05 / Actual Hollyford polygon · three bank gradients\nTanh · depth cap 10 m · assumed flat water level 0 m · 4 m grid',fontsize=16)
        fig.text(.5,.005,'Geometry-only experiment: all boundary edges are banks, including mouth/lake ends. No DEM, bed survey or flow simulation.\nLINZ river polygon: CC BY 4.0. Entire saved polygon shown; not clipped to the playable footprint.',ha='center',fontsize=9)
        fig.set_constrained_layout_pads(h_pad=.6,w_pad=.1)
        save(fig,out,'05-hollyford.png')
        metadata['examples']['hollyford']={'source':str(source.relative_to(ROOT)),'max_depth_m':10,'assumed_water_level_m':0,'resolution_m':4,'variants':stats}
    (out/'parameters.json').write_text(json.dumps(metadata,indent=2)+'\n')
    names=['01-profiles.png','02-lake-island.png','03-river.png','04-open-ends.png','05-hollyford.png']
    (out/'README.md').write_text('# Water-carver sample gallery\n\nGenerated directly by water-carver. Synthetic bed shapes, not surveyed bathymetry or game flow predictions.\n\n'+''.join(f'## {name[3:-4].replace("-"," ").title()}\n\n![{name}]({name})\n\n' for name in names if (out/name).exists()))
    print(out)
    print(json.dumps(metadata,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'output/water-carver-samples')
    render(parser.parse_args().output)
