"""Jamestown connected-centreline discovery and synthetic geometry trial.

This does not infer unattached streams or tributaries attached only to lakes.
All projected distances are metres. Source geometry remains read-only.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import numpy as np
import rasterio
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d,map_coordinates
import shapely
from shapely.geometry import LineString,MultiLineString,Point,Polygon,shape,mapping
from shapely.ops import transform,linemerge,unary_union,substring
from pyproj import Transformer
from nz_terrain import terrain
from nz_terrain.river.geometry import line_parts,transect,polygon_width
from widths import tapered_corridor
from joins import smooth_joins
from connector import connect
from smooth import fit_curve
from trial_mckenzie import extend,extend_into_polygons


def sample_dem(surface,aff,xy):
    xy=np.asarray(xy);c,r=(~aff)*(xy[:,0],xy[:,1])
    return map_coordinates(surface,[r-.5,c-.5],order=1,mode='constant',cval=np.nan)


def endpoint_extension(line,at_start,poly,distance=75.):
    station=0 if at_start else line.length
    near_station=min(10,line.length) if at_start else max(0,line.length-10)
    pt=np.array(line.interpolate(station).coords[0]);near=np.array(line.interpolate(near_station).coords[0])
    direction=(pt-near)/np.linalg.norm(pt-near)
    for angle in [0]+[sign*a for a in range(5,91,5) for sign in (1,-1)]:
        rad=np.deg2rad(angle);v=np.array([[np.cos(rad),-np.sin(rad)],[np.sin(rad),np.cos(rad)]])@direction
        if all(poly.covers(Point(pt+v*d)) for d in np.arange(2.5,35.1,2.5)):
            return pt+distance*v,angle
    raise ValueError('No 35 m interior sampling runway within 90 degrees')


def one_anchor(line,poly):
    """Free end first, polygon end last. No invented headwater relocation."""
    tip,turn=endpoint_extension(line,False,poly)
    route=LineString([*line.coords,tip]);samples=[]
    for offset in np.arange(2.5,75,5):
        s=line.length+offset;point=route.interpolate(s)
        if not poly.covers(point):continue
        # A side tributary's normal can run lengthwise down its receiving river.
        # Use the shortest complete local cross-section of that river polygon.
        options=[];centre=np.array(point.coords[0])
        for angle in np.deg2rad(np.arange(0,180,5)):
            direction=np.array([np.cos(angle),np.sin(angle)])
            section=LineString([centre-1000*direction,centre+1000*direction])
            width,status=polygon_width(section,point,poly)
            if width is None:continue
            piece=next(p for p in line_parts(section.intersection(poly)) if p.distance(point)<1e-6)
            options.append((width,piece))
        if not options:continue
        width,piece=min(options,key=lambda x:x[0])
        samples.append((offset,width,piece))
    if len(samples)<3:raise ValueError('Insufficient anchored width samples')
    selected=[s for s in samples if abs(s[0]-20)<=15]
    if len(selected)<3:raise ValueError('Insufficient samples around interior anchor')
    # Side tributaries keep their incoming direction. Midpoints of the receiving
    # channel's short cross-sections may coincide, so cannot define this tangent.
    position=np.array(route.interpolate(line.length+20).coords[0])
    slope=np.array(tip)-np.array(line.coords[-1]);slope/=np.linalg.norm(slope)
    length=line.length+20
    stations=np.linspace(0,length,max(3,int(length/2)+1))
    xy=np.array([route.interpolate(s).coords[0] for s in stations]);xy=gaussian_filter1d(xy,10/(stations[1]-stations[0]),axis=0,mode='nearest')
    inner=np.arange(40,line.length-20,40)
    start=np.array(line.coords[0]);direction=np.array(line.interpolate(min(10,line.length)).coords[0])-start;direction/=np.linalg.norm(direction)
    nodes=np.vstack([start,np.column_stack([np.interp(inner,stations,xy[:,k]) for k in (0,1)]),position])
    spline=CubicSpline(np.r_[0,inner,length],nodes,bc_type=((1,direction),(1,slope)))
    curve=LineString(spline(np.linspace(0,length,max(3,int(length)+1))))
    if not curve.is_simple:raise ValueError('One-anchor spline self-intersects')
    return curve,float(np.mean([s[1] for s in samples])),turn,samples


def cut_by_station(curve,widths,stations,keep):
    """Union cross-section strips only on accepted station intervals."""
    xy=np.array([curve.interpolate(s).coords[0] for s in stations])
    delta=np.array([np.array(curve.interpolate(min(curve.length,s+.5)).coords[0])-np.array(curve.interpolate(max(0,s-.5)).coords[0]) for s in stations])
    norm=np.linalg.norm(delta,axis=1);normal=np.column_stack([-delta[:,1],delta[:,0]])/norm[:,None]
    left=xy+normal*widths[:,None]/2;right=xy-normal*widths[:,None]/2
    return unary_union([tri for i in np.flatnonzero(keep[:-1]&keep[1:]) for tri in (Polygon([left[i],right[i],left[i+1]]),Polygon([right[i],right[i+1],left[i+1]]))])


def discover(root):
    loc=json.loads((root/'locations.json').read_text())['jamestown'];crs=terrain.projection(loc['longitude'],loc['latitude'])
    project=Transformer.from_crs(4326,crs,always_xy=True).transform
    def load(path):return [(f['properties'],transform(project,shape(f['geometry']))) for f in json.loads(path.read_text())['features']]
    polys={int(p['OBJECTID']):(p,g) for p,g in load(root/'river-configs/jamestown-polygons.geojson')}
    rivers=unary_union([g for p,g in polys.values()]);lakes=unary_union([g for p,g in load(root/'output/jamestown/metadata/lake-polygons.geojson')])
    features=load(root/'river-configs/jamestown-network.geojson')
    linz=[g for p,g in features if 'OpenStreetMap' not in p['source']]
    routes=[({'source':'LINZ merged degree-two lines'},g) for g in line_parts(linemerge(MultiLineString(linz)))]+[(p,g) for p,g in features if 'OpenStreetMap' in p['source']]
    endpoint_tree=shapely.STRtree([Point(g.coords[k]) for p,g in routes for k in (0,-1)])
    manifest=json.loads((root/'output/jamestown/metadata/manifest.json').read_text());cache=root/manifest['source_surface_cache']
    with rasterio.open(cache/'surface-world.tif') as ds:world=ds.read(1);world_aff=ds.transform
    footprint=Polygon([terrain.grid(loc,14336)*xy for xy in [(0,0),(4096,0),(4096,4096),(0,4096)]])
    candidates=[];skipped=[]
    for route_id,(properties,route) in enumerate(routes):
        for part_id,line in enumerate(line_parts(route.difference(rivers.union(lakes)))):
            if line.length<5 or not line.intersects(footprint):continue
            pts=[Point(line.coords[k]) for k in (0,-1)]
            ids=[min(polys,key=lambda i:polys[i][1].distance(p)) for p in pts]
            ids=[i if polys[i][1].distance(p)<=1 else None for i,p in zip(ids,pts)]
            item={'id':f'{route_id}-{part_id}','polygon_ids':ids,'route_properties':properties,'length_m':line.length}
            if all(i is None for i in ids):item['reason']='No river polygon anchor';skipped.append(item);continue
            if ids[0] is not None and ids[0]==ids[1]:item['reason']='Leaves and re-enters same river polygon; possible centreline drift';skipped.append(item);continue
            if all(i is not None for i in ids):item['kind']='polygon_connector'
            else:
                if ids[0] is not None:line=LineString(list(line.coords)[::-1]);pts.reverse();ids.reverse();item['polygon_ids']=ids
                end=pts[0];z=sample_dem(world,world_aff,[end.coords[0]])[0]
                if end.distance(lakes)<=1:item['kind']='water_link';item['water_end']='lake'
                elif np.isfinite(z) and z<0:item['kind']='water_link';item['water_end']='sea'
                elif len(endpoint_tree.query(end.buffer(1)))>1:item['kind']='branched';item['reason']='Free end is a mapped network junction'
                elif not footprint.buffer(-10).covers(end):item['kind']='boundary';item['reason']='Free end outside playable footprint; not an inferred headwater'
                else:item['kind']='tributary'
            candidates.append((item,line,route))
    return loc,crs,polys,rivers,lakes,footprint,cache,candidates,skipped


def generate(destination):
    root=terrain.BASE
    loc,crs,polys,rivers,lakes,footprint,cache,candidates,skipped=discover(root)
    with rasterio.open(cache/'surface-playable.tif') as ds:surface=ds.read(1);aff=ds.transform
    with rasterio.open(cache/'surface-world.tif') as ds:world=ds.read(1);waff=ds.transform
    inverse=Transformer.from_crs(crs,4326,always_xy=True).transform
    records=[];geometries=[];traces=[]
    for item,line,source_route in candidates:
        ids=item['polygon_ids'];kind=item['kind']
        if kind in ('branched','boundary'):
            skipped.append(item);continue
        try:
            if kind=='polygon_connector':
                pair=[polys[i][1] for i in ids]
                # Reuse the real interior route wherever both interior reaches exist.
                try:
                    connection=connect(source_route,*pair)
                    if connection.gap.hausdorff_distance(line)>1:raise ValueError('Different gap selected')
                    route=source_route
                except ValueError:
                    route=extend(line)
                    try:connection=connect(route,*pair)
                    except ValueError:
                        route,angles=extend_into_polygons(line,pair);item['extension_turn_degrees']=angles
                        connection=connect(route,*pair)
                curve,_,ends,_=fit_curve(route,connection)
                widths=connection.end_widths
                item['sample_widths_m']=list(widths)
            else:
                curve,width,angle,samples=one_anchor(line,polys[ids[1]][1]);item['sample_widths_m']=[width];item['extension_turn_degrees']=[angle];item['width_measurement']='Mean shortest local complete cross-section; 5 degree angular search'
                widths=(width,width) if kind=='water_link' else (.5,width)
            corridor,stations,widths_array,_=tapered_corridor(curve,widths,mode='linear')
            xy=np.array([curve.interpolate(s).coords[0] for s in stations]);z=sample_dem(surface,aff,xy)
            missing=~np.isfinite(z);z[missing]=sample_dem(world,waff,xy[missing])
            if not np.isfinite(z).all():raise ValueError('Route extends beyond available terrain')
            # Grade of a 20 m centred chord along the raw terrain profile; 10 m smoothing suppresses DEM noise.
            smoothed=gaussian_filter1d(z,10/(stations[1]-stations[0]),mode='nearest')
            lo=np.maximum(0,stations-10);hi=np.minimum(curve.length,stations+10)
            grade=np.degrees(np.arctan(np.abs((np.interp(hi,stations,smoothed)-np.interp(lo,stations,smoothed))/(hi-lo))))
            keep=(widths_array>=5)&(grade<15) if kind=='tributary' else np.ones(len(stations),bool)
            if kind=='tributary':corridor=cut_by_station(curve,widths_array,stations,keep)
            if corridor.is_empty:raise ValueError('No reach passes width and longitudinal slope thresholds')
            selected=[polys[i][1] for i in ids if i is not None]
            joined=smooth_joins(selected,corridor,20,60)
            patch=joined.synthetic_patch.difference(rivers.union(lakes)).intersection(footprint)
            if patch.is_empty:raise ValueError('No new area in playable footprint')
            item.update(status='generated',start_width_m=float(widths[0]),end_width_m=float(widths[1]),curve_length_m=curve.length,
                        passing_station_fraction=float(keep.mean()),max_grade_degrees=float(grade.max()),synthetic_area_m2=patch.area)
            # Preserve station evidence for raster-level threshold enforcement and review.
            traces.append({'id':item['id'],'kind':kind,'curve':curve,'stations':stations,'widths':widths_array,'grade':grade,'keep':keep})
            geometries.append((item,patch,curve));records.append(item)
            print(item['id'],kind,ids,'widths',np.round(widths,1),'area',round(patch.area),flush=True)
        except (ValueError,AssertionError) as exc:
            item['reason']=str(exc);item['status']='not_generated';skipped.append(item)
    destination.mkdir(parents=True,exist_ok=True)
    collection={'type':'FeatureCollection','features':[{'type':'Feature','properties':p,'geometry':mapping(transform(inverse,g))} for p,g,c in geometries]}
    (destination/'synthetic-connectors.geojson').write_text(json.dumps(collection)+'\n')
    (destination/'connector-traces.geojson').write_text(json.dumps({'type':'FeatureCollection','features':[{'type':'Feature','properties':{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in t.items() if k!='curve'},'geometry':mapping(transform(inverse,t['curve']))} for t in traces]})+'\n')
    report={'location':loc,'generated':records,'skipped':skipped,'parameters':{'headwater_width_m':.5,'minimum_tributary_carve_width_m':5,'maximum_tributary_grade_degrees':15,'grade_window_m':20,'grade_profile_smoothing_m':10,'width_interpolation':'linear','join_smoothing_m':20,'join_zone_m':60,'scope':'playable footprint; river polygon anchors only'}}
    (destination/'connection-report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('GENERATED',len(records),'SKIPPED',Counter(s.get('reason') for s in skipped),flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=terrain.BASE/'output/jamestown-river-lab/network-trial')
    args=parser.parse_args();generate(args.output)
