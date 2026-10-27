"""Addition-only Awarua seam trial; run after run_smooth.py."""
import argparse
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pyproj import Transformer
from shapely.geometry import shape,mapping
from shapely.ops import transform
from nz_terrain import terrain
from joins import smooth_joins


def parts(g):
    if g.is_empty:return
    if g.geom_type=='Polygon':yield g
    else:
        for part in getattr(g,'geoms',[]):yield from parts(part)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--from-id',type=int,default=2240)
    parser.add_argument('--to-id',type=int,default=4563)
    parser.add_argument('--radius',type=float,default=20.)
    parser.add_argument('--seam-reach',type=float,default=None)
    parser.add_argument('--curve-variant',choices=['smooth','variable-width','linear-width'],default='linear-width')
    args=parser.parse_args()
    if args.seam_reach is None: args.seam_reach=max(25.,3*args.radius)
    root=terrain.BASE;loc=json.loads((root/'locations.json').read_text())['jamestown']
    crs=terrain.projection(loc['longitude'],loc['latitude'])
    fwd=Transformer.from_crs(4326,crs,always_xy=True).transform
    inv=Transformer.from_crs(crs,4326,always_xy=True).transform
    pair=[transform(fwd,shape(f['geometry'])) for f in json.loads((root/'river-configs/jamestown-polygons.geojson').read_text())['features'] if int(f['properties']['OBJECTID']) in (args.from_id,args.to_id)]
    base=root/f'output/jamestown-river-lab/connectors/{args.from_id}-{args.to_id}'
    features=json.loads((base/args.curve_variant/'connector.geojson').read_text())['features']
    corridor=next(transform(fwd,shape(f['geometry'])) for f in features if f['properties']['kind']=='smooth-corridor')
    result=smooth_joins(pair,corridor,args.radius,args.seam_reach)
    out=base/(f'joins-{args.radius:g}m' if args.curve_variant=='smooth' else f'{args.curve_variant}-joins-{args.radius:g}m');out.mkdir(parents=True,exist_ok=True)
    report={'radius_m':args.radius,'seam_reach_m':args.seam_reach,'added_area_m2':result.additions.area,'removed_source_area_m2':result.original.difference(result.final).area,'removed_connector_area_m2':corridor.difference(result.final).area,'added_outside_seam_zone_m2':result.additions.difference(result.seam_zone).area,'reconstruction_error_m2':result.original.union(result.synthetic_patch).symmetric_difference(result.final).area,'valid':result.final.is_valid,'geometry_type':result.final.geom_type,'note':'Original polygons unchanged; all new water stored in synthetic patch. Closing only fills concave corners; no removals.'}
    assert report['valid']
    assert all(report[key]<1e-6 for key in ('removed_source_area_m2','removed_connector_area_m2','added_outside_seam_zone_m2','reconstruction_error_m2'))
    (out/'measurements.json').write_text(json.dumps(report,indent=2)+'\n')
    exported=[('synthetic-patch',result.synthetic_patch),('smoothing-additions',result.additions),('final-union',result.final)]
    (out/'connector.geojson').write_text(json.dumps({'type':'FeatureCollection','features':[{'type':'Feature','properties':{'kind':kind},'geometry':mapping(transform(inv,g))} for kind,g in exported]})+'\n')
    # Paired views share the exact scale and bounds at each end.
    fig,axes=plt.subplots(2,2,figsize=(13,10),layout='constrained')
    fit=json.loads((base/args.curve_variant/'measurements.json').read_text())
    limits=[(x-80,y-80,x+80,y+80) for x,y in fit['endpoint_xy']]
    for row,(ident,lim) in enumerate(zip(fit['endpoint_polygon_ids'],limits)):
        for col in range(2):
            ax=axes[row,col]
            for p in parts(result.original):ax.fill(*p.exterior.xy,color='#439fca',alpha=.6)
            for p in parts(corridor.difference(result.original)):ax.fill(*p.exterior.xy,color='#a67cce',alpha=.65)
            if col:
                for p in parts(result.additions):ax.fill(*p.exterior.xy,color='#f59e0b',zorder=5)
            boundary=result.before if col==0 else result.final
            for p in parts(boundary): ax.plot(*p.exterior.xy,color='#383442',lw=1.2,zorder=6)
            if col:
                for p in parts(result.before):ax.plot(*p.exterior.xy,color='#b44c31',ls=':',lw=1,zorder=7)
            ax.set(xlim=(lim[0],lim[2]),ylim=(lim[1],lim[3]),aspect='equal',title=f'{ident}: '+('before' if col==0 else f'addition-only smoothing ({args.radius:g} m)'),xlabel='Local easting (m)',ylabel='Local northing (m)');ax.grid(alpha=.18)
    fig.legend(handles=[Patch(fc='#439fca',label='Unchanged source polygons'),Patch(fc='#a67cce',label='Synthetic connector'),Patch(fc='#f59e0b',label='Added smoothing patches'),Line2D([],[],color='#b44c31',ls=':',label='Previous boundary')],loc='outside lower center',ncol=4)
    fig.suptitle(f'Awarua joins · {result.additions.area:.1f} m² added · no source or connector area removed',fontsize=15)
    fig.savefig(out/'comparison.png',dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,11),layout='constrained')
    for p in parts(result.original):ax.fill(*p.exterior.xy,color='#439fca',alpha=.6)
    for p in parts(corridor.difference(result.original)):ax.fill(*p.exterior.xy,color='#a67cce',alpha=.65)
    for p in parts(result.additions):ax.fill(*p.exterior.xy,color='#f59e0b',zorder=5)
    for p in parts(result.final):ax.plot(*p.exterior.xy,color='#383442',lw=.8,zorder=6)
    original_line=next(shape(f['geometry']) for f in features if f['properties']['kind']=='original-centreline')
    line=transform(fwd,original_line);ax.plot(*line.xy,ls=':',color='#b44c31',lw=1,label='Original centreline')
    a,b,c,d=corridor.bounds
    for ident,(x,y) in zip(fit['endpoint_polygon_ids'],fit['endpoint_xy']):ax.annotate(str(ident),(x,y),xytext=(25,15),textcoords='offset points',bbox={'fc':'white','ec':'none'},arrowprops={'arrowstyle':'-'})
    ax.set(xlim=(a-100,c+100),ylim=(b-100,d+100),aspect='equal',title=f'Awarua {args.from_id} ↔ {args.to_id} · {args.radius:g} m seam smoothing',xlabel='Local easting (m)',ylabel='Local northing (m)');ax.grid(alpha=.2)
    ax.legend(handles=[Patch(fc='#439fca',label='Unchanged source polygons'),Patch(fc='#a67cce',label='Synthetic connector'),Patch(fc='#f59e0b',label='Added smoothing'),Line2D([],[],ls=':',color='#b44c31',label='Original centreline')])
    fig.savefig(out/'overview.png',dpi=180);plt.close(fig)
    print(json.dumps(report,indent=2));print(out/'comparison.png')

if __name__=='__main__':main()
