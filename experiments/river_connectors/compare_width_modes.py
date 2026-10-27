"""Compare width interpolation on exactly the same saved Awarua curve."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from shapely.geometry import shape,mapping
from shapely.ops import transform
from pyproj import Transformer
from nz_terrain import terrain
from widths import tapered_corridor
from joins import smooth_joins
from run_joins import parts


def main():
    root=terrain.BASE;loc=json.loads((root/'locations.json').read_text())['jamestown']
    crs=terrain.projection(loc['longitude'],loc['latitude'])
    fwd=Transformer.from_crs(4326,crs,always_xy=True).transform
    inv=Transformer.from_crs(crs,4326,always_xy=True).transform
    base=root/'output/jamestown-river-lab/connectors/4563-1108'
    fit=json.loads((base/'variable-width/measurements.json').read_text())
    features=json.loads((base/'variable-width/connector.geojson').read_text())['features']
    curve=next(transform(fwd,shape(f['geometry'])) for f in features if f['properties']['kind']=='smooth-centreline')
    pair=[transform(fwd,shape(f['geometry'])) for f in json.loads((root/'river-configs/jamestown-polygons.geojson').read_text())['features'] if int(f['properties']['OBJECTID']) in (4563,1108)]
    out=base/'width-mode-comparison';out.mkdir(parents=True,exist_ok=True)
    results=[];reports=[]
    fig,ax=plt.subplots(figsize=(12,4),layout='constrained')
    labels=['200 m end transitions','Linear over entire reach']
    colours=['#70399b','#dc7625']
    for mode,label,colour in zip(('end-transitions','linear'),labels,colours):
        corridor,s,w,run=tapered_corridor(curve,fit['end_mean_widths_m'],mode=mode)
        result=smooth_joins(pair,corridor,20,60);results.append((corridor,result))
        ax.plot(s,w,label=label,color=colour,lw=2)
        report={'mode':mode,'start_width_m':float(w[0]),'end_width_m':float(w[-1]),'curve_length_m':curve.length,'corridor_area_m2':corridor.area,'seam_additions_m2':result.additions.area,'valid':result.final.is_valid,'source_removed_m2':result.original.difference(result.final).area,'reconstruction_error_m2':result.original.union(result.synthetic_patch).symmetric_difference(result.final).area}
        assert report['valid'] and report['source_removed_m2']<1e-6 and report['reconstruction_error_m2']<1e-6
        reports.append(report)
        (out/f'{mode}.geojson').write_text(json.dumps({'type':'FeatureCollection','features':[{'type':'Feature','properties':{'kind':kind,'width_mode':mode},'geometry':mapping(transform(inv,g))} for kind,g in [('synthetic-patch',result.synthetic_patch),('final-union',result.final),('smoothing-additions',result.additions)]]})+'\n')
    ax.set(xlabel='Distance along fitted centreline (m)',ylabel='Full channel width (m)',title='Same endpoints: 28.7 m → 18.2 m · same centreline');ax.grid(alpha=.2);ax.legend()
    fig.savefig(out/'width-profiles.png',dpi=180);plt.close(fig)
    def draw(ax,index,limits,title):
        corridor,result=results[index]
        for p in parts(result.original):ax.fill(*p.exterior.xy,color='#439fca',alpha=.65)
        for p in parts(corridor.difference(result.original)):ax.fill(*p.exterior.xy,color='#a67cce',alpha=.65)
        for p in parts(result.additions):ax.fill(*p.exterior.xy,color='#f59e0b',zorder=4)
        for p in parts(result.final):ax.plot(*p.exterior.xy,color='#383442',lw=.8,zorder=5)
        for p in parts(results[1-index][1].final):ax.plot(*p.exterior.xy,color='#bd3e56',ls=':',lw=1,zorder=6)
        ax.set(xlim=(limits[0],limits[2]),ylim=(limits[1],limits[3]),aspect='equal',title=title,xlabel='Local easting (m)',ylabel='Local northing (m)');ax.grid(alpha=.18)
    handles=[Patch(fc='#439fca',label='Source polygons'),Patch(fc='#a67cce',label='Connector'),Patch(fc='#f59e0b',label='20 m seam smoothing'),Line2D([],[],color='#bd3e56',ls=':',label='Other variant boundary')]
    fig,axes=plt.subplots(2,2,figsize=(13,10),layout='constrained')
    for row,((x,y),ident) in enumerate(zip(fit['endpoint_xy'],fit['endpoint_polygon_ids'])):
        for col in range(2):draw(axes[row,col],col,(x-80,y-80,x+80,y+80),f'{ident} · {labels[col]}')
    fig.legend(handles=handles,loc='outside lower center',ncol=4);fig.suptitle('Join comparison · identical scales and 20 m addition-only smoothing')
    fig.savefig(out/'joins.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,12),layout='constrained')
    a,b,c,d=curve.bounds
    for i,ax in enumerate(axes):draw(ax,i,(a-100,b-100,c+100,d+100),labels[i])
    fig.legend(handles=handles,loc='outside lower center',ncol=2)
    fig.savefig(out/'overview.png',dpi=180);plt.close(fig)
    (out/'measurements.json').write_text(json.dumps(reports,indent=2)+'\n')
    print(json.dumps(reports,indent=2))

if __name__=='__main__':main()
