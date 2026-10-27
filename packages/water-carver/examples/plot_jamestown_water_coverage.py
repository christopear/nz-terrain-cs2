"""Review Jamestown water polygons and uncovered centrelines; no terrain changes."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pyproj import Transformer
from shapely.geometry import shape, Polygon
from shapely.ops import transform, unary_union
from nz_terrain import terrain
from nz_terrain.river.geometry import line_parts


def main():
    root=terrain.BASE;out=root/'output/jamestown-river-lab/coverage';out.mkdir(parents=True,exist_ok=True)
    location=json.loads((root/'locations.json').read_text())['jamestown']
    project=Transformer.from_crs(4326,terrain.projection(location['longitude'],location['latitude']),always_xy=True).transform
    affine=terrain.grid(location,14336)
    footprint=Polygon([affine*p for p in [(0,0),(4096,0),(4096,4096),(0,4096)]])
    def load(path):
        return [(f['properties'],transform(project,shape(f['geometry'])).intersection(footprint)) for f in json.loads(path.read_text())['features']]
    rivers=[(p,g) for p,g in load(root/'river-configs/jamestown-polygons.geojson') if g.area>0]
    lakes=[(p,g) for p,g in load(root/'output/jamestown-polygon-rivers/metadata/lake-polygons.geojson') if g.area>0]
    water=unary_union([g for _,g in rivers+lakes])
    fig,ax=plt.subplots(figsize=(13,13),layout='constrained')
    for features,colour in [(lakes,'#b4dce8'),(rivers,'#398fbc')]:
        for p,g in features:
            for part in g.geoms if hasattr(g,'geoms') else [g]:
                if part.geom_type!='Polygon':continue
                x,y=part.exterior.xy;ax.fill(x,y,fc=colour,ec='#31768d',lw=.7,zorder=2)
                for hole in part.interiors:
                    x,y=hole.xy;ax.fill(x,y,fc='white',ec='#31768d',lw=.4,zorder=3)
    for p,g in load(root/'river-configs/jamestown-network.geojson'):
        for part in line_parts(g.difference(water)):
            x,y=part.xy;ax.plot(x,y,color='#c17635',lw=.85,zorder=4)
    names=['Hollyford','Awarua','McKenzie Creek','Hokuri Creek','Waiuna','McKerrow']
    for wanted in names:
        matches=[(p,g) for p,g in rivers+lakes if wanted.lower() in (p.get('name') or '').lower()]
        if not matches:continue
        p,g=max(matches,key=lambda item:item[1].area);pt=g.representative_point()
        ax.annotate(p['name'],(pt.x,pt.y),xytext=(15,14),textcoords='offset points',fontsize=10,
                    bbox=dict(fc='white',alpha=.9,ec='none'),arrowprops=dict(arrowstyle='-',color='.35'),zorder=6)
    x,y=footprint.exterior.xy;ax.plot(x,y,color='#333',lw=1.2)
    ax.set(aspect='equal',title='Jamestown / river and lake polygons + uncovered river centrelines\nRotated 14.336 km playable footprint · geometry review only',xlabel='Local easting (m)',ylabel='Local northing (m)')
    ax.legend(handles=[Patch(fc='#398fbc',label='River polygons'),Patch(fc='#b4dce8',label='Lake / lagoon polygons'),Line2D([],[],color='#c17635',label='Centrelines outside water polygons')],loc='lower left',framealpha=.95)
    ax.grid(alpha=.12)
    fig.savefig(out/'water-polygons-and-uncovered-centrelines.png',dpi=170);plt.close(fig)


if __name__=='__main__':main()
