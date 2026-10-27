"""Apply reviewed connected-centreline geometry to the existing Jamestown package."""
import argparse
from datetime import datetime,timezone
import hashlib,json,shutil
from pathlib import Path
import numpy as np
from PIL import Image
from affine import Affine
from pyproj import Transformer
import rasterio
from scipy.ndimage import gaussian_filter
import shapely
from shapely.geometry import shape,Point
from shapely.ops import transform,unary_union
from nz_terrain import terrain,hydrology
from nz_terrain.river.polygon_build import depth_field


def nearest_stations(points,curve,stations):
    """Exact projection onto sampled curve segments, avoiding vertex-only jumps."""
    xy=np.array([curve.interpolate(float(s)).coords[0] for s in stations])
    starts=xy[:-1];delta=xy[1:]-starts;norm=(delta*delta).sum(axis=1)
    segments=shapely.linestrings(np.stack([starts,xy[1:]],axis=1))
    tree=shapely.STRtree(segments)
    pairs=tree.query_nearest(points,all_matches=False)
    indices=np.empty(len(points),int);indices[pairs[0]]=pairs[1]
    coords=shapely.get_coordinates(points)
    frac=np.clip(((coords-starts[indices])*delta[indices]).sum(axis=1)/np.maximum(norm[indices],1e-12),0,1)
    return stations[indices]+frac*(stations[indices+1]-stations[indices])


def build(geometry_dir):
    root=terrain.BASE;slug='jamestown';label='jamestown-connected-rivers'
    base=root/'output/jamestown-polygon-rivers'
    loc=json.loads((root/'locations.json').read_text())[slug]
    baseline_manifest=json.loads((base/'metadata/baseline-manifest.json').read_text())
    if baseline_manifest['location']!=loc:raise ValueError('Baseline location mismatch')
    report=json.loads((geometry_dir/'connection-report.json').read_text())
    if report['location']!=loc:raise ValueError('Geometry location mismatch')
    settings=json.loads((root/'water-configs/river-preferred.json').read_text())
    crs=terrain.projection(loc['longitude'],loc['latitude']);project=Transformer.from_crs(4326,crs,always_xy=True).transform
    def read(path):return [(f['properties'],transform(project,shape(f['geometry']))) for f in json.loads(path.read_text())['features']]
    river_features=read(root/'river-configs/jamestown-polygons.geojson');original_water=unary_union([g for p,g in river_features])
    patches=read(geometry_dir/'synthetic-connectors.geojson')
    # Geographic round trips can invalidate sub-millimetre overlay slivers.
    # Repair synthetic patches only; retain authoritative polygon coordinates.
    patches=[(p,shapely.set_precision(shapely.make_valid(g),.001)) for p,g in patches]
    patch_union=unary_union([g for p,g in patches]);all_water=original_water.union(patch_union)
    trace_features=read(geometry_dir/'connector-traces.geojson');traces={p['id']:(p,g) for p,g in trace_features}
    opens=read(root/'river-configs/jamestown-hollyford-openings.geojson');opening_union=unary_union([g for p,g in opens])
    lake_open=next(g for p,g in opens if 'Lake' in p['name']);sea_open=next(g for p,g in opens if 'Sea' in p['name'])
    hollyford=unary_union([g for p,g in river_features if 'hollyford' in (p.get('name') or '').lower()])
    lake_features=json.loads((base/'metadata/lake-polygons.geojson').read_text())['features']
    lake_union=unary_union([transform(project,shape(f['geometry'])) for f in lake_features])
    # Recompute only the new channels and a local existing-river collar to remove old cutoff banks.
    collar=patch_union.buffer(60).intersection(original_water)
    components=list(all_water.geoms) if hasattr(all_water,'geoms') else [all_water]
    staging=root/'output'/f'.{label}-building';destination=root/'output'/label
    if staging.exists():shutil.rmtree(staging)
    staging.mkdir(parents=True)
    arrays={};changed={};audit={};eligible_counts={};component_report=[]
    for kind,width in terrain.WIDTHS.items():
        affine=terrain.grid(loc,width);spacing=width/4096
        old=np.array(Image.open(base/f'heightmaps/jamestown-polygon-rivers-{kind}.png'));result=old.copy()
        with rasterio.open(root/baseline_manifest['source_surface_cache']/f'surface-{kind}.tif') as ds:
            raw=ds.read(1)
            assert ds.shape==(4096,4096) and ds.transform.almost_equals(affine)
        allowed=np.zeros(old.shape,bool)
        window=hydrology.lake_window(collar,affine,4096)
        if window:
            r,c,mask=window;allowed[r:r+mask.shape[0],c:c+mask.shape[1]]|=mask
        for props,patch in patches:
            window=hydrology.lake_window(patch,affine,4096)
            if window is None:continue
            r,c,mask=window;h,w=mask.shape
            if props['kind']=='tributary' and mask.any():
                evidence,curve=traces[props['id']];rr,cc=np.where(mask);x,y=affine*(cc+c+.5,rr+r+.5)
                stations=np.asarray(evidence['stations'])
                along=nearest_stations(shapely.points(x,y),curve,stations)
                channel_width=np.interp(along,stations,evidence['widths']);grade=np.interp(along,stations,evidence['grade'])
                mask[rr,cc]=(channel_width>=5)&(grade<15)
            allowed[r:r+h,c:c+w]|=mask
            eligible_counts[f"{kind}:{props['id']}"]=int(mask.sum())
        # Lakes retain their established beds; no direct lake tributary excavation.
        lake_mask=hydrology.feature_mask(terrain,loc,kind,lake_features)
        allowed&=~lake_mask
        for number,water in enumerate(components):
            window=hydrology.lake_window(water,affine,4096)
            if window is None:continue
            r,c,mask=window;h,w=mask.shape;valid_domain=mask&allowed[r:r+h,c:c+w]
            if not valid_domain.any():continue
            local=affine*Affine.translation(c,r)
            # Opening edges are removed from banks even when the union changes their segmentation.
            banks=water.boundary.difference(opening_union.buffer(.05))
            print(kind,'component',number,'water pixels',int(mask.sum()),'eligible',int(valid_domain.sum()),flush=True)
            depth=depth_field(water,banks,local,mask.shape,crs,settings)
            values=np.maximum(raw[r:r+h,c:c+w],0)
            weights=gaussian_filter(mask.astype(float),30/spacing)
            reference=gaussian_filter(np.where(mask,values,0),30/spacing)/np.maximum(weights,1e-12)
            rr,cc=np.indices(mask.shape);x,y=local*(cc+.5,rr+.5)
            holly=shapely.contains_xy(hollyford,x,y)
            if holly.any():
                points=shapely.points(x[holly],y[holly])
                dl=shapely.distance(points,lake_open);ds=shapely.distance(points,sea_open)
                reference[holly]=ds/(dl+ds)
            valid=valid_domain&np.isfinite(depth)
            encoded=terrain.encode(np.where(valid,reference-depth,0))
            view=result[r:r+h,c:c+w];view[valid]=np.minimum(view[valid],encoded[valid])
            component_report.append({'grid':kind,'component':number,'eligible_pixels':int(valid.sum()),'maximum_synthetic_depth_m':float(depth[valid].max())})
        assert np.array_equal(result[~allowed],old[~allowed])
        assert np.all(result<=old)
        arrays[kind]=result;changed[kind]=result!=old
        audit[kind]={'allowed_pixels':int(allowed.sum()),'changed_pixels_before_nesting':int(changed[kind].sum()),'outside_allowed_changed_pixels':int((changed[kind]&~allowed).sum()),'max_additional_cut_m':float((old.astype(float)-result).max()*terrain.STEP),'lake_pixels_changed':int((changed[kind]&lake_mask).sum())}
        print(kind,audit[kind],flush=True)
    affected=changed['playable'].reshape(1024,4,1024,4).any(axis=(1,3))
    means=arrays['playable'].astype(float).reshape(1024,4,1024,4).mean(axis=(1,3));centre=arrays['world'][1536:2560,1536:2560]
    centre[affected]=np.rint(means[affected]).astype('uint16')
    error=float(np.abs(centre[affected]-means[affected]).max()*terrain.STEP) if affected.any() else 0
    assert error<=terrain.STEP/2+1e-10
    for kind,array in arrays.items():
        path=staging/f'{label}-{kind}.png';Image.fromarray(array).save(path)
        check=np.array(Image.open(path));assert check.dtype==np.uint16 and check.shape==(4096,4096);np.testing.assert_array_equal(array,check)
        old=np.array(Image.open(base/f'heightmaps/jamestown-polygon-rivers-{kind}.png'));audit[kind]['changed_pixels_exported']=int((array!=old).sum())
        # Include actual river masks in preview colouring; previews do not imply water sources exist.
        river_mask=np.zeros(array.shape,bool)
        win=hydrology.lake_window(all_water,terrain.grid(loc,terrain.WIDTHS[kind]),4096)
        if win:
            r,c,m=win;river_mask[r:r+m.shape[0],c:c+m.shape[1]]=m
        terrain.preview(array,loc,kind,staging/f'{kind}-preview.jpg',lake_mask=hydrology.feature_mask(terrain,loc,kind,lake_features)|river_mask)
    for name in ['synthetic-connectors.geojson','connector-traces.geojson','connection-report.json']:
        shutil.copy(geometry_dir/name,staging/name)
    for name in ['lake-polygons.geojson','water-report.json','baseline-manifest.json','jamestown-polygons.geojson','jamestown-hollyford-openings.geojson','SOURCES.md']:
        if (base/'metadata'/name).exists():shutil.copy(base/'metadata'/name,staging/name)
    shutil.copy(root/'water-configs/river-preferred.json',staging/'river-preferred.json')
    report_out={'generated_utc':datetime.now(timezone.utc).isoformat(),'location':loc,'baseline':'jamestown-polygon-rivers','baseline_sha256':{k:hashlib.sha256((base/f'heightmaps/jamestown-polygon-rivers-{k}.png').read_bytes()).hexdigest() for k in terrain.WIDTHS},'geometry_sha256':hashlib.sha256((geometry_dir/'synthetic-connectors.geojson').read_bytes()).hexdigest(),'audit':audit,'nested_grid_max_error_m':error,'eligible_pixels_by_connector':eligible_counts,'components':component_report,'settings':settings,'limitations':['Synthetic bed stages use a 30 m smoothed raw DEM; no hydraulic simulation or monotonic-flow guarantee.','Tributary grade means absolute longitudinal slope of a 20 m chord on a 10 m smoothed raw-terrain profile.','Mapped river polygons and between-polygon connectors retain their prior unrestricted slope policy.','Existing deep valleys are not specially detected; lower absolute bed target wins.','Width comes from receiving river geometry, not a surveyed tributary width.','Source water is preserved; addition-only smoothing can fill nearby land holes.','Only playable-footprint synthetic additions are included; the world map centre is synchronised.','Water sources must be placed manually in-game.']}
    report_out['source_snapshots_sha256']={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'river-configs/jamestown-network.geojson',root/'river-configs/jamestown-polygons.geojson',root/'output/jamestown/metadata/lake-polygons.geojson']}
    report_out['scope']=report.get('scope','Connected-river network trial')
    report_out['included_connectors']=len(report['generated'])
    (staging/'connected-river-report.json').write_text(json.dumps(report_out,indent=2)+'\n')
    instructions=(base/'IMPORT.txt').read_text().replace('jamestown-polygon-rivers',label)
    instructions+='\nCONNECTED RIVER TRIAL\nUse BOTH heightmaps from this package. River-polygon connectors use linear endpoint widths. See connection-report.json for the selected scope and exclusions.\nThis is an experimental geometry and terrain package, not hydraulically validated. Place water sources manually.\n'
    instructions+='\nSELECTED SCOPE\n'+report.get('scope','Connected-river network trial')+'\n'
    (staging/'IMPORT.txt').write_text(instructions)
    archive=terrain.package_outputs(staging,label)
    if destination.exists():shutil.rmtree(destination)
    staging.replace(destination)
    print('COMPLETE',destination/archive.name,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--geometry-dir',type=Path,default=terrain.BASE/'output/jamestown-river-lab/reviewed-connectors')
    args=parser.parse_args();build(args.geometry_dir)
