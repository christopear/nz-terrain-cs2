"""Compare addition-only seam radii, preserving the earlier individual trial."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pyproj import Transformer
from shapely.geometry import shape,mapping,box
from shapely.ops import transform
from nz_terrain import terrain
from joins import smooth_joins
from run_joins import parts


def main():
    root=terrain.BASE;loc=json.loads((root/'locations.json').read_text())['jamestown']
    crs=terrain.projection(loc['longitude'],loc['latitude'])
    fwd=Transformer.from_crs(4326,crs,always_xy=True).transform
    inv=Transformer.from_crs(crs,4326,always_xy=True).transform
    pair=[transform(fwd,shape(f['geometry'])) for f in json.loads((root/'river-configs/jamestown-polygons.geojson').read_text())['features'] if int(f['properties']['OBJECTID']) in (2240,4563)]
    base=root/'output/jamestown-river-lab/connectors/2240-4563'
    features=json.loads((base/'smooth/connector.geojson').read_text())['features']
    corridor=next(transform(fwd,shape(f['geometry'])) for f in features if f['properties']['kind']=='smooth-corridor')
    out=base/'join-radius-sweep';out.mkdir(parents=True,exist_ok=True)
    results=[];reports=[]
    for radius in (5,10,25,50):
        reach=max(25,3*radius)
        result=smooth_joins(pair,corridor,radius,reach);results.append(result)
        report={'radius_m':radius,'seam_reach_m':reach,'added_area_m2':result.additions.area,'removed_source_area_m2':result.original.difference(result.final).area,'removed_connector_area_m2':corridor.difference(result.final).area,'reconstruction_error_m2':result.original.union(result.synthetic_patch).symmetric_difference(result.final).area,'added_outside_zone_m2':result.additions.difference(result.seam_zone).area,'valid':result.final.is_valid,'added_area_within_1cm_of_zone_edge_m2':result.additions.intersection(result.seam_zone.boundary.buffer(.01)).area}
        assert report['valid']
        assert all(report[k]<1e-6 for k in ('removed_source_area_m2','removed_connector_area_m2','reconstruction_error_m2','added_outside_zone_m2'))
        reports.append(report)
        (out/f'joins-{radius}m.geojson').write_text(json.dumps({'type':'FeatureCollection','features':[{'type':'Feature','properties':{'kind':kind,'radius_m':radius,'seam_reach_m':reach},'geometry':mapping(transform(inv,g))} for kind,g in [('synthetic-patch',result.synthetic_patch),('smoothing-additions',result.additions),('final-union',result.final)]]})+'\n')
    (out/'measurements.json').write_text(json.dumps(reports,indent=2)+'\n')
    limits=[(2720,4410,2900,4590),(2490,4640,2670,4820)]
    def draw(ax,result,lim,title):
        for p in parts(result.original):ax.fill(*p.exterior.xy,color='#439fca',alpha=.6)
        for p in parts(corridor.difference(result.original)):ax.fill(*p.exterior.xy,color='#a67cce',alpha=.65)
        for p in parts(result.additions):ax.fill(*p.exterior.xy,color='#f59e0b',zorder=5)
        for p in parts(result.final):ax.plot(*p.exterior.xy,color='#383442',lw=1,zorder=6)
        for p in parts(result.before):ax.plot(*p.exterior.xy,color='#b44c31',ls=':',lw=.8,zorder=7)
        ax.set(xlim=(lim[0],lim[2]),ylim=(lim[1],lim[3]),aspect='equal',title=title,xlabel='Local easting (m)',ylabel='Local northing (m)');ax.grid(alpha=.18);ax.tick_params(labelsize=8)
    handles=[Patch(fc='#439fca',label='Unchanged source'),Patch(fc='#a67cce',label='Connector'),Patch(fc='#f59e0b',label='Added fill'),Line2D([],[],color='#b44c31',ls=':',label='Before smoothing')]
    fig,axes=plt.subplots(2,4,figsize=(20,10),layout='constrained')
    for row,(ident,lim) in enumerate(zip((2240,4563),limits)):
        for col,(result,report) in enumerate(zip(results,reports)):
            local_area=result.additions.intersection(box(*lim)).area
            draw(axes[row,col],result,lim,f"{ident} · radius {report['radius_m']} m\n{local_area:.1f} m² added in view")
    fig.legend(handles=handles,loc='outside lower center',ncol=4)
    fig.suptitle('Awarua · addition-only smoothing · equal 180 × 180 m views\nRadii 5 / 10 / 25 / 50 m; local seam zones 25 / 30 / 75 / 150 m',fontsize=15)
    fig.savefig(out/'comparison.png',dpi=180);plt.close(fig)
    # A whole-connector view exposes any broader effects outside the detail panels.
    fig,axes=plt.subplots(1,4,figsize=(20,7),layout='constrained')
    for ax,result,report in zip(axes,results,reports):draw(ax,result,(2400,4380,2920,4850),f"Radius {report['radius_m']} m · {report['added_area_m2']:.1f} m² total")
    fig.legend(handles=handles,loc='outside lower center',ncol=4)
    fig.savefig(out/'overview.png',dpi=170);plt.close(fig)
    print(json.dumps(reports,indent=2))

if __name__=='__main__':main()
