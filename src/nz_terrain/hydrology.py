"""Synthetic water beds from independent absolute elevations, never additive cuts."""
from datetime import datetime, timezone
import json
import hashlib
from pathlib import Path
from urllib.parse import urlencode
import numpy as np
from affine import Affine
from pyproj import Transformer
from rasterio.features import rasterize
from scipy.ndimage import distance_transform_edt, gaussian_filter
from shapely.geometry import shape, mapping, Polygon
from shapely.ops import transform

LAKE_SERVICE='https://services.arcgis.com/xdsHIIxuCWByZiCB/ArcGIS/rest/services/NZ_Lake_Polygons_Topo_1_50k/FeatureServer/0'
LAGOON_SERVICE='https://services6.arcgis.com/ZVM1rEuVZjtC1Wwk/ArcGIS/rest/services/Hydrography_Data/FeatureServer/18'


def combine_targets(base, targets):
    """Replace terrain where a target owns it; choose lower bed in overlaps.

    Targets are (absolute elevation array, boolean domain mask), not depths to
    subtract. The unchanged land is not a competing target inside water masks.
    """
    result=base.copy();owned=np.zeros(base.shape,bool)
    for target,mask in targets:
        if target.shape!=base.shape or mask.shape!=base.shape:
            raise ValueError('Water target grid mismatch')
        if not np.isfinite(target[mask]).all():raise ValueError('Nonfinite bed target')
        overlap=mask & owned;fresh=mask & ~owned
        result[fresh]=target[fresh]
        result[overlap]=np.minimum(result[overlap],target[overlap])
        owned|=mask
    return result


def fetch_lakes(terrain,slug,location,kind='lake'):
    service=LAKE_SERVICE if kind=='lake' else LAGOON_SERVICE
    bbox=terrain.bounds(terrain.outline(location,57344))
    key=hashlib.sha256(json.dumps(bbox).encode()).hexdigest()[:12]
    directory=terrain.BASE/'data'/'water'/slug/key;directory.mkdir(parents=True,exist_ok=True)
    path=directory/f'{kind}s.geojson'
    if path.exists():return json.loads(path.read_text())
    params={'f':'geojson','where':'1=1','geometry':','.join(map(str,bbox)),
            'geometryType':'esriGeometryEnvelope','inSR':4326,'spatialRel':'esriSpatialRelIntersects',
            'outFields':'*','outSR':4326,'returnGeometry':'true','resultRecordCount':1000,
            'orderByFields':'OBJECTID' if kind=='lake' else 'FID'}
    features=[];offset=0
    while True:
        params['resultOffset']=offset
        page=terrain.fetch(service+'/query?'+urlencode(params))
        if 'error' in page:raise ValueError(f'Lake service error: {page["error"]}')
        if 'features' not in page:raise ValueError('Lake service returned no feature collection')
        batch=page['features'];features.extend(batch);offset+=len(batch)
        if not page.get('exceededTransferLimit') and not page.get('properties',{}).get('exceededTransferLimit') and len(batch)<1000:break
        if not batch:raise ValueError('Lake service pagination stalled')
    ids=[f['properties']['t50_fid'] for f in features]
    if len(ids)!=len(set(ids)):raise ValueError('Lake service repeated a page')
    data={'type':'FeatureCollection','features':features,'source':service,
          'retrieved_utc':datetime.now(timezone.utc).isoformat(),
          'attribution':'LINZ Topo50 lake polygons, CC BY 4.0; public ArcGIS distribution.'}
    temporary=path.with_suffix('.tmp');temporary.write_text(json.dumps(data));temporary.replace(path)
    return data


def lake_window(geometry,affine,size,padding=0):
    """Rasterise a local lake window, retaining polygon islands and rotated grids."""
    inv=~affine
    pixel=transform(lambda x,y,z=None:inv*(x,y),geometry)
    x0,y0,x1,y1=pixel.bounds
    c0=max(-padding,int(np.floor(x0))-2);c1=min(size+padding,int(np.ceil(x1))+2)
    r0=max(-padding,int(np.floor(y0))-2);r1=min(size+padding,int(np.ceil(y1))+2)
    if c1<=c0 or r1<=r0:return None
    mask=rasterize([(mapping(geometry),1)],out_shape=(r1-r0,c1-c0),
                   transform=affine*Affine.translation(c0,r0),fill=0,dtype='uint8').astype(bool)
    return r0,c0,mask


def interior_window(window,size):
    r,c,mask=window
    r0=max(0,r);c0=max(0,c);r1=min(size,r+mask.shape[0]);c1=min(size,c+mask.shape[1])
    return (slice(r0,r1),slice(c0,c1)),(slice(r0-r,r1-r),slice(c0-c,c1-c))


def basin_depth(mask,spacing,maximum,dropoff,smoothing,shore_transition=0.,shore_gradient=.01):
    if not all(np.isfinite(v) for v in [spacing,maximum,dropoff,smoothing,shore_transition,shore_gradient]) or spacing<=0 or maximum<=0 or dropoff<=0 or smoothing<0 or shore_transition<0 or shore_gradient<0:
        raise ValueError('Invalid lake depth/dropoff/smoothing/shore parameters')
    # Half-cell correction puts zero depth at the raster shoreline, rather than
    # one complete pixel inland. Islands are holes and become shore boundaries.
    distance=np.maximum(distance_transform_edt(mask,sampling=spacing)-spacing/2,0).astype('float32')
    if smoothing:
        distance=np.minimum(distance,gaussian_filter(distance,smoothing/spacing,mode='constant'))
    if shore_transition:
        from water_carver.river import profile
        return profile(distance,shore_transition,maximum,shore_gradient,maximum/dropoff)
    return maximum*np.tanh(distance/dropoff)


def water_settings(terrain,slug,options=None):
    defaults=json.loads((terrain.BASE/'water-configs'/'defaults.json').read_text())
    local=terrain.BASE/'water-configs'/f'{slug}.json'
    if local.exists():
        override=json.loads(local.read_text());overrides={**defaults['overrides'],**override.get('overrides',{})}
        defaults.update(override);defaults['overrides']=overrides
    defaults.update(options or {})
    for kind in ['sea','lake','lagoon']:
        for key in [f'{kind}_depth_m',f'{kind}_gradient']:
            if not np.isfinite(defaults[key]) or defaults[key]<=0:
                raise ValueError(f'{key} must be finite and positive')
    for key in ['lagoon_distance_smoothing_m','lagoon_shore_transition_m','lagoon_shore_gradient']:
        if key in defaults and (not np.isfinite(defaults[key]) or defaults[key]<0):
            raise ValueError(f'{key} must be finite and nonnegative')
    return defaults


def prepare_lakes(terrain,slug,location,surface,options=None):
    defaults=water_settings(terrain,slug,options)
    data=fetch_lakes(terrain,slug,location)
    lagoons=fetch_lakes(terrain,slug,location,'lagoon')
    features={str(f['properties']['t50_fid']):f for f in data['features']}
    for feature in lagoons['features']:
        feature['properties']['water_kind']='lagoon'
        features[str(feature['properties']['t50_fid'])]=feature
    supplement_path=terrain.BASE/'water-configs'/f'{slug}-polygons.geojson'
    supplement=None
    if supplement_path.exists():
        supplement=json.loads(supplement_path.read_text())
        for feature in supplement['features']:
            features[str(feature['properties']['t50_fid'])]=feature
    crs=terrain.projection(location['longitude'],location['latitude'])
    project=Transformer.from_crs(4326,crs,always_xy=True).transform
    world_affine=terrain.grid(location,57344)
    footprint=Polygon([world_affine*p for p in [(0,0),(4096,0),(4096,4096),(0,4096)]])
    prepared=[];skipped=[]
    for feature in features.values():
        p=feature['properties'];identifier=str(p['t50_fid']);override=defaults['overrides'].get(identifier,{})
        kind=override.get('kind',p.get('water_kind','lagoon' if 'lagoon' in (p.get('name') or '').lower() else 'lake'))
        geometry=transform(project,shape(feature['geometry']))
        if not geometry.is_valid:raise ValueError(f'Invalid lake polygon {identifier}; repair explicitly')
        if not geometry.intersects(footprint):continue
        if override.get('skip') or p.get('perennial') in {'dry','seasonal'} or geometry.area<defaults['minimum_lake_area_m2']:
            skipped.append({'id':identifier,'name':p.get('name'),'reason':'Excluded by config, area, or seasonal/dry classification'});continue
        window=lake_window(geometry,world_affine,4096)
        if window is None or not window[2].any():
            # A small lake may miss every coarse world-grid pixel centre while
            # still covering playable pixels. Keep its polygon and estimate its
            # level from the nearest source sample rather than dropping it.
            point=geometry.intersection(footprint).representative_point()
            c,r=(~world_affine)*(point.x,point.y)
            values=np.array([surface['world'][np.clip(int(r),0,4095),np.clip(int(c),0,4095)]])
        else:
            slices,sub=interior_window(window,4096);values=surface['world'][slices][window[2][sub]]
        observed=float(np.median(values));attribute=p.get('elevation')
        if 'water_level_m' in override:level=float(override['water_level_m']);level_source='explicit override'
        elif attribute is not None and np.isfinite(attribute) and abs(float(attribute)-observed)<=5:
            level=float(attribute);level_source='LINZ water-surface elevation'
        else:level=observed;level_source='median uncarved DEM inside lake (estimate)'
        maximum=float(override.get('depth_m',defaults[f'{kind}_depth_m']))
        dropoff=float(override.get('dropoff_m',maximum/defaults[f'{kind}_gradient']))
        if not np.isfinite(level) or level-maximum<terrain.HEIGHT_MIN or level>terrain.HEIGHT_MIN+terrain.HEIGHT_SPAN:
            raise ValueError(f'Lake {identifier} outside heightmap range')
        prepared.append({'geometry':geometry,'feature':feature,'report':{
            'id':identifier,'name':p.get('name') or f'Unnamed {kind} {identifier}','kind':kind,'area_m2':geometry.area,
            'polygon_source':p.get('source',data['source']),'polygon_source_url':p.get('source_url'),
            'water_level_m':level,'water_level_source':level_source,'dem_median_m':observed,
            'dem_p10_p90_m':np.percentile(values,[10,90]).tolist(),
            'maximum_depth_parameter_m':maximum,'depth_source':override.get('depth_source','default synthetic depth'),
            'dropoff_m':dropoff,
            'distance_smoothing_m':float(override.get('distance_smoothing_m',defaults.get(f'{kind}_distance_smoothing_m',defaults['distance_smoothing_m']))),
            'shore_transition_m':float(override.get('shore_transition_m',defaults.get(f'{kind}_shore_transition_m',0.))),
            'shore_gradient':float(override.get('shore_gradient',defaults.get(f'{kind}_shore_gradient',.01))),
            'grids':{}}})
    masks={};targets={}
    for name,width in terrain.WIDTHS.items():
        affine=terrain.grid(location,width);spacing=width/4096
        targets[name]=surface[name].copy();masks[name]=np.zeros((4096,4096),bool)
        for lake in prepared:
            record=lake['report'];padding=int(np.ceil(6*record['dropoff_m']/spacing))+2
            # First check intersection with this grid before allocating padding.
            test=lake_window(lake['geometry'],affine,4096)
            if test is None or not test[2].any():continue
            window=lake_window(lake['geometry'],affine,4096,padding)
            mask=window[2];depth=basin_depth(mask,spacing,record['maximum_depth_parameter_m'],record['dropoff_m'],record['distance_smoothing_m'],record['shore_transition_m'],record['shore_gradient'])
            dst,src=interior_window(window,4096);inside=mask[src]
            bed=record['water_level_m']-depth[src]
            view=targets[name][dst];owned=masks[name][dst];overlap=inside & owned
            view[inside & ~owned]=bed[inside & ~owned];view[overlap]=np.minimum(view[overlap],bed[overlap]);owned|=inside
            record['grids'][name]={'pixels':int(inside.sum()),'deepest_synthetic_depth_m':float(depth[src][inside].max()),
                                  'lowest_bed_m':float(bed[inside].min())}
        print(f'{slug} {name}: {int(masks[name].sum())} lake pixels',flush=True)
    report={'method':'Independent target = water level minus synthetic depth from smoothed shore distance. A configured shoreline transition smoothly ramps the initial slope to max_depth/dropoff before the tanh cap. Lower absolute target wins overlaps; no additive depth.',
            'defaults':defaults,'source':data['source'],'lagoon_source':lagoons['source'],'retrieved_utc':data['retrieved_utc'],
            'lakes':[x['report'] for x in prepared],'skipped':skipped,
            'limitations':['All basin geometry is synthetic, even where depth parameters have a published source.',
                           'Small/narrow lakes may not approach their maximum depth parameter.',
                           'Water levels with no credible attribute use the original DEM median and need review.',
                           'Mapped perennial lake and lagoon polygons are carved; separate pond layers and unmapped water are not automatically included.',
                           'No features does not prove no water: source coverage can be incomplete, especially offshore islands.',
                           'Heightmaps do not place water sources, control lake levels, or calibrate flow.']}
    geojson={'type':'FeatureCollection','features':[x['feature'] for x in prepared],
             'attribution':data['attribution']}
    if supplement is not None:
        report['supplement_sources']=supplement.get('source_urls',[])
        geojson['attribution']+=' Supplementary polygons: © OpenStreetMap contributors, ODbL 1.0; https://www.openstreetmap.org/copyright. Adapted polygon database supplied in this package.'
    return targets,masks,report,geojson


def feature_mask(terrain,location,name,features):
    project=Transformer.from_crs(4326,terrain.projection(location['longitude'],location['latitude']),always_xy=True).transform
    shapes=[(mapping(transform(project,shape(f['geometry']))),1) for f in features]
    if not shapes:return np.zeros((4096,4096),bool)
    return rasterize(shapes,out_shape=(4096,4096),transform=terrain.grid(location,terrain.WIDTHS[name]),dtype='uint8').astype(bool)


def plot_water_profiles(terrain,location,raw,final,report,features,out):
    """Before/after transects through deepest generated cells; no simulated water."""
    import os
    os.environ.setdefault('MPLCONFIGDIR','/tmp/jamestown-matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    available=[r for r in report['lakes'] if 'world' in r['grids']]
    if not available:return
    selected=sorted(available,key=lambda r:(not any(n in r['name'] for n in ['McKerrow','Waiuna']),-r['area_m2']))[:3]
    by_id={str(f['properties']['t50_fid']):f for f in features['features']}
    project=Transformer.from_crs(4326,terrain.projection(location['longitude'],location['latitude']),always_xy=True).transform
    affine=terrain.grid(location,57344)
    fig,axes=plt.subplots(len(selected),1,figsize=(11,3.5*len(selected)),squeeze=False,layout='constrained')
    for ax,record in zip(axes.flat,selected):
        geometry=transform(project,shape(by_id[record['id']]['geometry']))
        window=lake_window(geometry,affine,4096);dst,sub=interior_window(window,4096)
        mask=window[2][sub];z=final['world'][dst].astype('float32')*terrain.STEP+terrain.HEIGHT_MIN
        deepest=np.argmin(np.where(mask,z,np.inf));r,c=np.unravel_index(deepest,mask.shape)
        columns=np.flatnonzero(mask[r]);left=max(0,columns.min()-10);right=min(mask.shape[1],columns.max()+11)
        old=raw['world'][dst][r,left:right].astype('float32')*terrain.STEP+terrain.HEIGHT_MIN
        new=z[r,left:right];distance=(np.arange(left,right)-c)*14
        ax.plot(distance,old,label='Original DEM surface',color='.5')
        ax.plot(distance,new,label='Synthetic bed',color='tab:blue')
        ax.axhline(record['water_level_m'],color='teal',ls='--',label='Intended water level (manual editor setup)')
        ax.set(title=record['name'],xlabel='Across deepest world-grid row (m)',ylabel='Elevation (m)',
               ylim=(record['water_level_m']-record['maximum_depth_parameter_m']*1.1,record['water_level_m']+10))
        ax.grid(alpha=.2)
    axes.flat[0].legend(fontsize=8);fig.savefig(out/'water-profiles.png',dpi=140);plt.close(fig)
