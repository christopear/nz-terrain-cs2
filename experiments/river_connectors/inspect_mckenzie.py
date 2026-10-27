"""Inspect McKenzie polygon gaps and divergent centreline sources, without carving."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from pyproj import Transformer
from shapely.geometry import shape,box
from shapely.ops import transform,unary_union
from nz_terrain import terrain
from nz_terrain.river.geometry import line_parts
from run_joins import parts


def main():
    root=terrain.BASE;loc=json.loads((root/'locations.json').read_text())['jamestown']
    project=Transformer.from_crs(4326,terrain.projection(loc['longitude'],loc['latitude']),always_xy=True).transform
    def read(path):return [(f['properties'],transform(project,shape(f['geometry']))) for f in json.loads((root/path).read_text())['features']]
    polygons=[(p,g) for p,g in read('river-configs/jamestown-polygons.geojson') if 'mckenzie' in p.get('name','').lower()]
    union=unary_union([g for p,g in polygons]);window=box(-50,1700,2550,3150)
    lines=[(p,g) for p,g in read('river-configs/jamestown-network.geojson') if g.intersects(window)]
    report=[]
    for p,g in lines:
        nearby=g.intersection(window)
        report.append({'name':p.get('name'),'id':p.get('OBJECTID',p.get('river_section_id')),'source':p.get('source'),'length_in_view_m':nearby.length,'inside_polygons_m':nearby.intersection(union).length,'polygon_distances_m':{str(q['OBJECTID']):g.distance(poly) for q,poly in polygons}})
    out=root/'output/jamestown-river-lab/connectors/mckenzie-review';out.mkdir(parents=True,exist_ok=True)
    (out/'source-review.json').write_text(json.dumps(report,indent=2)+'\n')
    fig=plt.figure(figsize=(16,10),layout='constrained');grid=fig.add_gridspec(2,3,height_ratios=[1.2,1])
    axes=[fig.add_subplot(grid[0,:])]+[fig.add_subplot(grid[1,i]) for i in range(3)]
    limits=[(-50,1700,2550,3150),(580,2550,830,2800),(1430,2430,1680,2680),(2220,1850,2470,2100)]
    titles=['McKenzie Creek · full source geometry','180 ↔ 2238 · 250 m view','2238 ↔ 2239 · 250 m view','2239 ↔ 4562 · 250 m view']
    for ax,lim,title in zip(axes,limits,titles):
        view=box(*lim)
        for p,g in polygons:
            for part in parts(g):ax.fill(*part.exterior.xy,fc='#439fca',ec='#236583',alpha=.6)
            if g.intersects(view):
                pt=g.intersection(view).representative_point();ax.annotate(str(p['OBJECTID']),(pt.x,pt.y),xytext=(0,16),textcoords='offset points',bbox={'fc':'white','ec':'none','alpha':.8},fontsize=9)
        for p,g in lines:
            named='mckenzie' in p.get('name','').lower()
            for part in line_parts(g.intersection(view)):
                ax.plot(*part.xy,color='#228b45' if named else '#888888',lw=2 if named else 1,ls='-' if named else '--',zorder=4)
        ax.set(xlim=(lim[0],lim[2]),ylim=(lim[1],lim[3]),aspect='equal',title=title,xlabel='Local easting (m)',ylabel='Local northing (m)');ax.grid(alpha=.18)
    fig.legend(handles=[Patch(fc='#439fca',label='LINZ river polygons'),Line2D([],[],color='#228b45',label='Named McKenzie LINZ lines'),Line2D([],[],color='#888888',ls='--',label='Other LINZ lines')],loc='outside lower center',ncol=3)
    fig.suptitle('McKenzie Creek · source comparison before any new connectors',fontsize=16)
    fig.savefig(out/'source-comparison.png',dpi=180);plt.close(fig)
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
