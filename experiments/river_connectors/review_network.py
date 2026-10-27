"""Review geometry and exported connected-river heightmaps with independent audits."""
import json,zipfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import shapely
from PIL import Image
from pyproj import Transformer
from shapely.geometry import shape,Polygon
from shapely.ops import transform
from nz_terrain import terrain,hydrology
from run_joins import parts
from trial_mckenzie import fill_polygon


def main():
    root=terrain.BASE;out=root/'output/jamestown-connected-rivers';meta=out/'metadata';previews=out/'previews'
    report=json.loads((meta/'connection-report.json').read_text());loc=report['location'];crs=terrain.projection(loc['longitude'],loc['latitude'])
    fwd=Transformer.from_crs(4326,crs,always_xy=True).transform
    def read(path):return [(f['properties'],transform(fwd,shape(f['geometry']))) for f in json.loads(path.read_text())['features']]
    polys=read(meta/'jamestown-polygons.geojson');lakes=read(meta/'lake-polygons.geojson');patches=read(meta/'synthetic-connectors.geojson');traces=read(meta/'connector-traces.geojson')
    affine=terrain.grid(loc,14336);footprint=Polygon([affine*p for p in [(0,0),(4096,0),(4096,4096),(0,4096)]])
    colours={'polygon_connector':'#8b54bc','tributary':'#d57b28','water_link':'#dc4a89'}
    def draw(ax,limits,title):
        for p,g in lakes:
            for part in parts(g):fill_polygon(ax,part,'#abd8ea',.7)
        for p,g in polys:
            for part in parts(g):fill_polygon(ax,part,'#439fca',.7)
        for p,g in traces:
            if p['kind']=='tributary':ax.plot(*g.xy,color='#888888',ls=':',lw=.8,zorder=1)
        for p,g in patches:
            for part in parts(g):fill_polygon(ax,part,colours[p['kind']],.85)
        ax.plot(*footprint.exterior.xy,color='.4',ls='--',lw=.7)
        ax.set(xlim=(limits[0],limits[2]),ylim=(limits[1],limits[3]),aspect='equal',title=title,xlabel='Local easting (m)',ylabel='Local northing (m)');ax.grid(alpha=.15)
    handles=[Patch(fc='#439fca',label='River polygons'),Patch(fc='#abd8ea',label='Lakes / lagoons'),Patch(fc=colours['polygon_connector'],label='Two-polygon connectors'),*[Patch(fc=colours[k],label=k.replace('_',' ').title()) for k in ('tributary','water_link') if any(p['kind']==k for p,g in patches)]]
    fig,ax=plt.subplots(figsize=(12,12),layout='constrained');draw(ax,footprint.bounds,'Jamestown · reviewed Awarua and McKenzie connectors')
    fig.legend(handles=handles,loc='outside lower center',ncol=2)
    fig.savefig(previews/'river-network-overview.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(15,13),layout='constrained')
    windows=[((1100,3400,4600,7000),'Big Bay / Awarua'),((-100,600,2800,3200),'McKenzie Creek'),((-7600,-2200,-4400,2300),'Lower Hollyford connections'),((-1600,-7600,1500,-4500),'Hokuri Creek connections')]
    for ax,(bounds,title) in zip(axes.flat,windows):draw(ax,bounds,title)
    fig.legend(handles=handles,loc='outside lower center',ncol=3)
    fig.suptitle('Reviewed polygon connectors only · linear widths and 20 m join smoothing',fontsize=15)
    fig.savefig(previews/'river-network-details.png',dpi=180);plt.close(fig)
    # Actual export difference, sharing a scale for both panels.
    old=np.array(Image.open(root/'output/jamestown-polygon-rivers/heightmaps/jamestown-polygon-rivers-playable.png'))
    new=np.array(Image.open(out/'heightmaps/jamestown-connected-rivers-playable.png'))
    delta=(old.astype(float)-new)*terrain.STEP
    fig,ax=plt.subplots(figsize=(10,10),layout='constrained')
    z=old.astype(float)*terrain.STEP+terrain.HEIGHT_MIN
    ax.imshow(z[::4,::4],cmap='Greys',vmin=0,vmax=1500,extent=(0,14.336,14.336,0),alpha=.4)
    cuts=delta.reshape(1024,4,1024,4).max(axis=(1,3))
    im=ax.imshow(np.ma.masked_less_equal(cuts,0),cmap='magma',vmin=0,vmax=max(10,delta.max()),extent=(0,14.336,14.336,0))
    ax.set(title='Actual additional excavation in exported playable heightmap',xlabel='Map x (km)',ylabel='Map y (km)');fig.colorbar(im,ax=ax,label='Additional lowering (m)',shrink=.7)
    fig.savefig(previews/'additional-excavation.png',dpi=170);plt.close(fig)
    checks={};arrays={}
    for kind in terrain.WIDTHS:
        path=out/f'heightmaps/jamestown-connected-rivers-{kind}.png';array=np.array(Image.open(path));assert array.dtype==np.uint16 and array.shape==(4096,4096)
        original=np.array(Image.open(root/f'output/jamestown-polygon-rivers/heightmaps/jamestown-polygon-rivers-{kind}.png'))
        arrays[kind]=array;checks[kind]={'dtype':str(array.dtype),'shape':list(array.shape),'changed_pixels':int((array!=original).sum()),'raised_pixels':int((array>original).sum())}
        if kind=='playable':assert checks[kind]['raised_pixels']==0
    affected=(new!=old).reshape(1024,4,1024,4).any(axis=(1,3));means=new.astype(float).reshape(1024,4,1024,4).mean(axis=(1,3))
    error=float(np.abs(arrays['world'][1536:2560,1536:2560][affected]-means[affected]).max()*terrain.STEP)
    assert error<=terrain.STEP/2+1e-10
    checks['nesting_max_error_m']=error
    # Independent check against vector source boundaries and exact curve projection.
    rr,cc=np.where(new!=old);x,y=affine*(cc+.5,rr+.5);points=shapely.points(x,y)
    source_union=shapely.union_all([g for p,g in polys])
    covered=shapely.covers(source_union,points)
    trace_by_id={p['id']:(p,g) for p,g in traces}
    for p,g in patches:
        idx=np.flatnonzero(~covered & shapely.covers(shapely.make_valid(g).buffer(.002),points))
        if not len(idx):continue
        if p['kind']=='tributary':
            evidence,curve=trace_by_id[p['id']]
            station=shapely.line_locate_point(curve,points[idx])
            widths=np.interp(station,evidence['stations'],evidence['widths'])
            grades=np.interp(station,evidence['stations'],evidence['grade'])
            ok=(widths>=5-1e-5)&(grades<15+1e-5)
            covered[idx[ok]]=True
        else:covered[idx]=True
    checks['changed_pixels_outside_source_or_eligible_synthetic_water']=int((~covered).sum())
    assert covered.all(),checks
    lake_union=shapely.union_all([g for p,g in lakes])
    checks['playable_lake_pixels_changed_vector_audit']=int(shapely.covers(lake_union,points).sum())
    assert checks['playable_lake_pixels_changed_vector_audit']==0
    (meta/'export-validation.json').write_text(json.dumps(checks,indent=2)+'\n')
    archive=terrain.package_outputs(out,'jamestown-connected-rivers')
    with zipfile.ZipFile(archive) as z:assert z.testzip() is None
    print(json.dumps(checks,indent=2));print(archive)

if __name__=='__main__':main()
