"""Big Bay geometry review: complete centrelines and mapped water polygons."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pyproj import Transformer
from shapely.geometry import shape, box
from shapely.ops import transform
from nz_terrain import terrain
from nz_terrain.river.geometry import line_parts


def main():
    root=terrain.BASE;out=root/'output/jamestown-river-lab/coverage'
    loc=json.loads((root/'locations.json').read_text())['jamestown']
    project=Transformer.from_crs(4326,terrain.projection(loc['longitude'],loc['latitude']),always_xy=True).transform
    def load(path):
        return [(f['properties'],transform(project,shape(f['geometry']))) for f in json.loads(path.read_text())['features']]
    rivers=load(root/'river-configs/jamestown-polygons.geojson')
    lakes=load(root/'output/jamestown-polygon-rivers/metadata/lake-polygons.geojson')
    lines=load(root/'river-configs/jamestown-network.geojson')
    bounds=(1200,3500,5200,7400)
    fig,axes=plt.subplots(1,2,figsize=(16,9),layout='constrained')
    for ax,limits,title in zip(axes,[bounds,(1800,3900,3400,5700)],['Big Bay / Awarua River and Waiuna Lagoon','Closer view / intermittent Awarua polygons']):
        window=box(*limits)
        for features,colour in [(lakes,'#b7deeb'),(rivers,'#439fca')]:
            for p,g in features:
                if not g.intersects(window):continue
                clipped=g.intersection(window)
                for part in clipped.geoms if hasattr(clipped,'geoms') else [clipped]:
                    if part.geom_type!='Polygon':continue
                    x,y=part.exterior.xy;ax.fill(x,y,fc=colour,ec='#236583',lw=1,zorder=2)
                    for hole in part.interiors:
                        x,y=hole.xy;ax.fill(x,y,fc='white',ec='#236583',lw=.6,zorder=3)
        for p,g in lines:
            for part in line_parts(g.intersection(window)):
                x,y=part.xy;ax.plot(x,y,color='#c75a20',lw=1.4,zorder=4)
        # Label separate mapped Awarua polygon features, rather than implying continuity.
        awarua=[(p,g) for p,g in rivers if 'awarua' in (p.get('name') or '').lower()]
        for p,g in awarua:
            if not g.intersects(window):continue
            pt=g.intersection(window).representative_point()
            ax.annotate(f'Awarua polygon {p["OBJECTID"]}',(pt.x,pt.y),xytext=(-75,15),textcoords='offset points',fontsize=9,bbox=dict(fc='white',alpha=.9,ec='none'),arrowprops=dict(arrowstyle='-',color='.35'),zorder=6)
        if ax is axes[0]:
            ax.text(4050,5450,'Waiuna Lagoon',ha='center',fontsize=11,bbox=dict(fc='white',alpha=.8,ec='none'),zorder=6)
        ax.set(xlim=(limits[0],limits[2]),ylim=(limits[1],limits[3]),aspect='equal',title=title,xlabel='Local easting (m)',ylabel='Local northing (m)');ax.grid(alpha=.18)
    axes[0].legend(handles=[Patch(fc='#439fca',label='River polygons'),Patch(fc='#b7deeb',label='Lake / lagoon polygons'),Line2D([],[],color='#c75a20',label='Full centrelines (including inside polygons)')],loc='upper left',fontsize=9)
    fig.suptitle('Big Bay / full mapped centrelines over water polygons\nOriginal source geometry, before width/slope filtering · no inferred connections or terrain changes',fontsize=15)
    fig.savefig(out/'big-bay-full-centrelines.png',dpi=180);plt.close(fig)


if __name__=='__main__':main()
