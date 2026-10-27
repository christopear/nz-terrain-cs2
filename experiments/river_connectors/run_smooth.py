"""Compare the original connector with a midpoint-aligned smooth curve."""
import argparse
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import shape,mapping,Point
from shapely.ops import transform,substring
from pyproj import Transformer
from nz_terrain import terrain
from connector import connect
from smooth import fit_curve
from widths import tapered_corridor


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-id",type=int,default=2240)
    parser.add_argument("--to-id",type=int,default=4563)
    parser.add_argument("--width-transition",type=float,default=None,help="Enable endpoint widths; transition length in metres")
    parser.add_argument("--width-mode",choices=["linear","end-transitions","constant"],default="linear")
    args=parser.parse_args()
    ids=(args.from_id,args.to_id)
    root=terrain.BASE
    loc=json.loads((root/'locations.json').read_text())['jamestown']
    crs=terrain.projection(loc['longitude'],loc['latitude'])
    fwd=Transformer.from_crs(4326,crs,always_xy=True).transform
    inv=Transformer.from_crs(crs,4326,always_xy=True).transform
    polygons={int(f['properties']['OBJECTID']):transform(fwd,shape(f['geometry'])) for f in json.loads((root/'river-configs/jamestown-polygons.geojson').read_text())['features']}
    pair=[polygons[i] for i in ids]
    routes=[transform(fwd,shape(f['geometry'])) for f in json.loads((root/'river-configs/jamestown-network.geojson').read_text())['features']]
    matches=[r for r in routes if all(r.intersects(p) for p in pair)]
    if len(matches)!=1: raise ValueError('Expected one shared route')
    route=matches[0];base=connect(route,*pair)
    params=dict(anchor_depth=20.,fit_radius=15.,guide_spacing=40.,guide_smoothing=10.,output_spacing=1.)
    curve,corridor,ends,guides=fit_curve(route,base,**params)
    variant={'constant':'smooth','end-transitions':'variable-width','linear':'linear-width'}[args.width_mode]
    out=root/f'output/jamestown-river-lab/connectors/{ids[0]}-{ids[1]}/{variant}';out.mkdir(parents=True,exist_ok=True)
    if args.width_mode != "constant":
        corridor,distances,widths,run=tapered_corridor(curve,base.end_widths,args.width_transition or 200.,mode=args.width_mode)
        fig,ax=plt.subplots(figsize=(12,3),layout='constrained')
        ax.plot(distances,widths,color='#70399b');ax.set(xlabel='Distance along fitted curve (m)',ylabel='Full channel width (m)',title=f'Width interpolation: {args.width_mode} · {run:g} m');ax.grid(alpha=.2)
        fig.savefig(out/'width-profile.png',dpi=180);plt.close(fig)
    source=substring(route,ends[0][0],ends[1][0])
    offsets=[route.distance(curve.interpolate(s)) for s in np.linspace(0,curve.length,1000)]
    end_ids=[ids[min(range(2),key=lambda i:pair[i].distance(Point(e[1])))] for e in ends]
    report={'polygon_ids':ids,'endpoint_polygon_ids':end_ids,'endpoint_xy':[e[1].tolist() for e in ends],'gap_length_m':base.gap.length,'end_mean_widths_m':base.end_widths,'parameters_m':params,'width_m':base.width,'width_transition_m':args.width_transition,'width_mode':args.width_mode,'curve_length_m':curve.length,'matching_source_length_m':source.length,'max_distance_to_source_m':max(offsets),'mean_distance_to_source_m':float(np.mean(offsets)),'valid_corridor':bool(corridor.is_valid),'simple_curve':bool(curve.is_simple),'intersects_both_polygons':all(corridor.intersects(p) for p in pair),'note':'Width mode recorded above; mean width is the interior target. Anchors are 20m inside source polygons; endpoint directions fitted from nearby transect midpoints. Not bank blending or a carving change.'}
    (out/'measurements.json').write_text(json.dumps(report,indent=2)+'\n')
    features=[{'type':'Feature','properties':{'kind':kind},'geometry':mapping(transform(inv,g))} for kind,g in [('smooth-centreline',curve),('smooth-corridor',corridor),('original-centreline',source)]]
    (out/'connector.geojson').write_text(json.dumps({'type':'FeatureCollection','features':features})+'\n')
    fig,axes=plt.subplots(2,2,figsize=(14,12),layout='constrained')
    a,b,c,d=base.gap.bounds
    limits=[(a-90,b-90,c+90,d+90)]*2
    limits +=[(e[1][0]-90,e[1][1]-90,e[1][0]+90,e[1][1]+90) for e in ends]
    titles=['Original: buffered mapped centreline','Trial: aligned endpoints + smooth guided curve',f'Join detail: polygon {end_ids[0]}',f'Join detail: polygon {end_ids[1]}']
    for n,(ax,lim,title) in enumerate(zip(axes.flat,limits,titles)):
        for p in pair:
            for part in p.geoms if hasattr(p,'geoms') else [p]: ax.fill(*part.exterior.xy,color='#439fca',alpha=.65)
        bed=base.corridor if n==0 else corridor
        ax.fill(*bed.exterior.xy,color='#a67cce',alpha=.4)
        ax.plot(*route.xy,color='#c75a20',lw=1.3,ls='--',label='Mapped line')
        if n: ax.plot(*curve.xy,color='#70399b',lw=2,label='Smooth line')
        for e in ends:
            ax.scatter(*e[3].T,s=12,color='#16624b',zorder=5)
            ax.scatter(*e[1],s=55,marker='x',color='black',zorder=6)
        ax.set(xlim=(lim[0],lim[2]),ylim=(lim[1],lim[3]),aspect='equal',title=title,xlabel='Local easting (m)',ylabel='Local northing (m)');ax.grid(alpha=.2)
        ax.legend(loc='best',fontsize=8)
    fig.suptitle(f'Awarua {ids[0]} ↔ {ids[1]} · width mode: {report["width_mode"]}\nBlue: source polygons · purple: connector · green: fitted transect midpoints · ×: interior anchors',fontsize=13)
    fig.savefig(out/'comparison.png',dpi=170);plt.close(fig)
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
