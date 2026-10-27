"""Treat each McKenzie gap line independently, extending its end tangents."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch,PathPatch
from matplotlib.path import Path
from shapely.geometry.polygon import orient
from shapely.geometry import shape,mapping,LineString,Point,box
from shapely.ops import transform,unary_union
from pyproj import Transformer
from nz_terrain import terrain
from connector import connect
from smooth import fit_curve
from widths import tapered_corridor
from joins import smooth_joins
from run_joins import parts


def fill_polygon(ax,poly,colour,alpha=.6):
    poly=orient(poly,sign=1.)
    vertices=[];codes=[]
    for ring in [poly.exterior,*poly.interiors]:
        xy=list(ring.coords);vertices.extend(xy)
        codes.extend([Path.MOVETO]+[Path.LINETO]*(len(xy)-2)+[Path.CLOSEPOLY])
    ax.add_patch(PathPatch(Path(vertices,codes),facecolor=colour,edgecolor='none',alpha=alpha))


def extend(line,distance=75.,tangent_run=10.):
    points=np.array(line.coords)
    start=np.array(line.interpolate(0).coords[0]);near=np.array(line.interpolate(min(tangent_run,line.length)).coords[0])
    end=np.array(line.interpolate(line.length).coords[0]);before=np.array(line.interpolate(max(0,line.length-tangent_run)).coords[0])
    a=start-near;b=end-before
    if min(np.linalg.norm(a),np.linalg.norm(b))<1e-9:raise ValueError('Degenerate endpoint tangent')
    return LineString(np.vstack([start+distance*a/np.linalg.norm(a),points,end+distance*b/np.linalg.norm(b)]))



def extend_into_polygons(line, polygons, distance=75., tangent_run=10.):
    """Use straight extensions where possible; minimally turn misses inward.

    Search in 5 degree increments, requiring 35 m of interior sampling runway.
    The observed gap line is untouched; only invented interior extensions turn.
    """
    xy=np.array(line.coords); tips=[];turns=[]
    for station,near_station in ((0,min(tangent_run,line.length)),(line.length,max(0,line.length-tangent_run))):
        pt=np.array(line.interpolate(station).coords[0]);near=np.array(line.interpolate(near_station).coords[0])
        direction=(pt-near)/np.linalg.norm(pt-near)
        poly=min(polygons,key=lambda p:p.distance(Point(pt)))
        selected=None
        for angle in [0]+[sign*a for a in range(5,91,5) for sign in (1,-1)]:
            rad=np.deg2rad(angle);rot=np.array([[np.cos(rad),-np.sin(rad)],[np.sin(rad),np.cos(rad)]])
            v=rot@direction
            if all(poly.covers(Point(pt+v*d)) for d in np.arange(2.5,35.1,2.5)):
                selected=(pt+distance*v,angle);break
        if selected is None:raise ValueError('No inward extension within 90 degrees')
        tips.append(selected[0]);turns.append(selected[1])
    return LineString(np.vstack([tips[0],xy,tips[1]])),turns


def main():
    root=terrain.BASE;loc=json.loads((root/'locations.json').read_text())['jamestown']
    crs=terrain.projection(loc['longitude'],loc['latitude'])
    fwd=Transformer.from_crs(4326,crs,always_xy=True).transform
    inv=Transformer.from_crs(crs,4326,always_xy=True).transform
    polygons={int(f['properties']['OBJECTID']):transform(fwd,shape(f['geometry'])) for f in json.loads((root/'river-configs/jamestown-polygons.geojson').read_text())['features'] if int(f['properties']['OBJECTID']) in (180,2238,2239,4562)}
    lines={int(f['properties']['OBJECTID']):transform(fwd,shape(f['geometry'])) for f in json.loads((root/'river-configs/jamestown-network.geojson').read_text())['features'] if f['properties'].get('OBJECTID') in (380913,677434,677435,677436,677437)}
    groups=[((180,2238),(380913,677434)),((2238,2239),(677435,)),((2239,4562),(677436,677437))]
    out=root/'output/jamestown-river-lab/connectors/mckenzie-separate-lines';out.mkdir(parents=True,exist_ok=True)
    reports=[];results={};features=[]
    colours=['#a67cce','#e29d49']
    for pair,ids in groups:
        for ident in ids:
            line=lines[ident];route=extend(line)
            report={'line_id':ident,'polygon_ids':pair,'extension_m':75,'tangent_run_m':10,'original_line_length_m':line.length,'extended_interior_lengths_m':[route.intersection(polygons[i]).length for i in pair]}
            features.append(('extended-line',ident,route))
            try:
                try:
                    connection=connect(route,*[polygons[i] for i in pair])
                except ValueError as exc:
                    report['straight_extension_failure']=str(exc)
                    route,turns=extend_into_polygons(line,[polygons[i] for i in pair])
                    report['interior_extension_turn_degrees']=turns
                    features.append(('adjusted-extended-line',ident,route))
                    connection=connect(route,*[polygons[i] for i in pair])
                curve,_,ends,_=fit_curve(route,connection)
                corridor,s,w,_=tapered_corridor(curve,connection.end_widths)
                result=smooth_joins([polygons[i] for i in pair],corridor,20,60)
                assert result.final.is_valid and result.original.difference(result.final).area<1e-6
                assert result.original.union(result.synthetic_patch).symmetric_difference(result.final).area<1e-6
                report.update(status='generated',end_widths_route_order_m=connection.end_widths,sample_counts=[sum(s['end']==i for s in connection.samples) for i in (0,1)],smoothing_added_m2=result.additions.area)
                results[ident]=(result,curve,connection)
                features.extend([('synthetic-patch',ident,result.synthetic_patch),('smooth-line',ident,curve)])
                features.extend(('transect',ident,s['geometry']) for s in connection.samples)
            except (ValueError,AssertionError) as exc:
                report.update(status='not_generated',reason=str(exc))
            reports.append(report)
    original=unary_union(list(polygons.values()))
    final=original.union(unary_union([r[0].synthetic_patch for r in results.values()]))
    features.append(('combined-final',None,final))
    (out/'geometry.geojson').write_text(json.dumps({'type':'FeatureCollection','features':[{'type':'Feature','properties':{'kind':kind,'line_id':ident},'geometry':mapping(transform(inv,g))} for kind,ident,g in features]})+'\n')
    (out/'measurements.json').write_text(json.dumps(reports,indent=2)+'\n')
    fig,axes=plt.subplots(2,3,figsize=(16,11),layout='constrained')
    limits=[(555,2500,855,2800),(1395,2400,1695,2700),(2190,1800,2490,2100)]
    for col,((pair,ids),lim) in enumerate(zip(groups,limits)):
        for row in (0,1):
            ax=axes[row,col]
            for ident in pair:
                for p in parts(polygons[ident]):fill_polygon(ax,p,'#439fca',.55)
            for j,ident in enumerate(ids):
                if ident in results:
                    result,curve,connection=results[ident]
                    if row==0:
                        for p in parts(result.synthetic_patch):fill_polygon(ax,p,colours[j],.55)
                        for sample in connection.samples:ax.plot(*sample['geometry'].xy,color='#236b40',alpha=.55,lw=.7)
                        ax.plot(*curve.xy,color=colours[j],lw=1.4)
                ax.plot(*lines[ident].xy,color='#333333',ls=':',lw=1.2)
                if row==0:
                    pt=lines[ident].interpolate(.5,normalized=True);ax.annotate(str(ident),(pt.x,pt.y),xytext=(5,8),textcoords='offset points',fontsize=8,bbox={'fc':'white','ec':'none','alpha':.8})
            if row==1:
                patch=unary_union([results[i][0].synthetic_patch for i in ids if i in results])
                for p in parts(patch):fill_polygon(ax,p,'#a67cce',.65)
                joined=unary_union([polygons[i] for i in pair]).union(patch)
                for p in parts(joined):
                    for ring in [p.exterior,*p.interiors]:ax.plot(*ring.xy,color='#383442',lw=1)
            ax.set(xlim=(lim[0],lim[2]),ylim=(lim[1],lim[3]),aspect='equal',title=f'{pair[0]} ↔ {pair[1]} · '+('separate channels + transects' if row==0 else 'combined water footprint'),xlabel='Local easting (m)',ylabel='Local northing (m)');ax.grid(alpha=.2)
    fig.legend(handles=[Patch(fc='#439fca',label='Source polygons'),Patch(fc=colours[0],label='Channel 1'),Patch(fc=colours[1],label='Channel 2')],loc='outside lower center',ncol=3)
    fig.suptitle('McKenzie · independent extended lines · linear widths · 20 m addition-only joins\nAll panels show 300 × 300 m; dotted lines are the source routes',fontsize=14)
    fig.savefig(out/'comparison.png',dpi=180);plt.close(fig)
    assert final.is_valid and original.difference(final).area<1e-6
    print(json.dumps(reports,indent=2))

if __name__=='__main__':main()
