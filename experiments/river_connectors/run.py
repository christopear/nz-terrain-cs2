"""Run with uv run python experiments/river_connectors/run.py."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from pyproj import Transformer
from shapely.geometry import shape, mapping
from shapely.ops import transform
from nz_terrain import terrain
from connector import connect


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--from-id',type=int,default=2240)
    parser.add_argument('--to-id',type=int,default=4563)
    parser.add_argument('--sample-run',type=float,default=75)
    parser.add_argument('--sample-step',type=float,default=5)
    parser.add_argument('--tangent-window',type=float,default=10)
    args=parser.parse_args()
    root=terrain.BASE
    loc=json.loads((root/'locations.json').read_text())['jamestown']
    crs=terrain.projection(loc['longitude'],loc['latitude'])
    forward=Transformer.from_crs(4326,crs,always_xy=True).transform
    inverse=Transformer.from_crs(crs,4326,always_xy=True).transform
    polygons={int(f['properties']['OBJECTID']):transform(forward,shape(f['geometry'])) for f in json.loads((root/'river-configs/jamestown-polygons.geojson').read_text())['features']}
    pair=[polygons[i] for i in (args.from_id,args.to_id)]
    routes=[]
    for f in json.loads((root/'river-configs/jamestown-network.geojson').read_text())['features']:
        g=transform(forward,shape(f['geometry']))
        if all(g.intersects(p) for p in pair): routes.append((f['properties'],g))
    if len(routes)!=1: raise ValueError(f'Expected one shared source route; found {len(routes)}')
    props,route=routes[0]
    result=connect(route,*pair,args.sample_run,args.sample_step,args.tangent_window)
    out=root/f'output/jamestown-river-lab/connectors/{args.from_id}-{args.to_id}'
    out.mkdir(parents=True,exist_ok=True)
    report={'polygon_ids':[args.from_id,args.to_id],'source_route':props,'parameters':vars(args),'gap_length_m':result.gap.length,'end_mean_widths_route_order_m':result.end_widths,'connector_width_m':result.width,'samples':[{k:v for k,v in s.items() if k!='geometry'} for s in result.samples],'note':'Synthetic constant-width proposal; original centreline retained; no terrain edits. Widths are plan-view transects, not depth measurements.'}
    (out/'measurements.json').write_text(json.dumps(report,indent=2)+'\n')
    features=[]
    for kind,g in [('connector',result.corridor),('gap-centreline',result.gap)]+[('width-sample',s['geometry']) for s in result.samples]:
        features.append({'type':'Feature','properties':{'kind':kind},'geometry':mapping(transform(inverse,g))})
    (out/'connector.geojson').write_text(json.dumps({'type':'FeatureCollection','features':features})+'\n')
    fig,ax=plt.subplots(figsize=(11,9),layout='constrained')
    for ident,poly in zip((args.from_id,args.to_id),pair):
        for p in poly.geoms if hasattr(poly,'geoms') else [poly]:
            ax.fill(*p.exterior.xy,color='#439fca',alpha=.7,zorder=2)
        pt=poly.intersection(result.gap.buffer(130)).representative_point()
        ax.annotate(str(ident),(pt.x,pt.y),xytext=(30,25),textcoords='offset points',arrowprops={'arrowstyle':'-'},bbox={'fc':'white','ec':'none'},zorder=6)
    ax.fill(*result.corridor.exterior.xy,fc='#a67cce',ec='#70399b',alpha=.45,zorder=3)
    ax.plot(*route.xy,color='#c75a20',lw=1.5,zorder=4)
    for s in result.samples: ax.plot(*s['geometry'].xy,color='#16624b',lw=1,alpha=.8,zorder=5)
    a,b,c,d=result.gap.bounds
    ax.set(xlim=(a-130,c+130),ylim=(b-130,d+130),aspect='equal',xlabel='Local easting (m)',ylabel='Local northing (m)',title=f'Awarua: {args.from_id} ↔ {args.to_id}\nMean endpoint widths {result.end_widths[0]:.1f} / {result.end_widths[1]:.1f} m → {result.width:.1f} m connector · gap {result.gap.length:.0f} m')
    ax.legend(handles=[Patch(fc='#439fca',label='Original river polygons'),Patch(fc='#a67cce',label='Proposed constant-width connector'),Line2D([],[],color='#c75a20',label='Original mapped centreline'),Line2D([],[],color='#16624b',label='Perpendicular width samples')],loc='lower left')
    ax.grid(alpha=.2)
    fig.savefig(out/'connector-preview.png',dpi=180)
    print(json.dumps({k:v for k,v in report.items() if k!='samples'},indent=2))
    print(out/'connector-preview.png')

if __name__=='__main__': main()
