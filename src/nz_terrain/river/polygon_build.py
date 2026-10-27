"""Build a preset polygon-river trial over a matching sea/lake baseline.

No centreline carving. Reference levels outside the lower Hollyford are DEM
estimates; this command produces an experimental package for in-game review.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

from affine import Affine
import numpy as np
from PIL import Image
from pyproj import Transformer
import rasterio
from scipy.ndimage import gaussian_filter
import shapely
from shapely.geometry import LineString, MultiLineString, Point, shape
from shapely.ops import transform
from water_carver.river import profile, smooth_creases
from .. import hydrology, terrain


def retained_banks(water, openings):
    """Remove only complete reviewed exterior edges; retain all island shores."""
    parts=list(water.geoms) if water.geom_type=='MultiPolygon' else [water]
    edges=[];removed=[]
    for part in parts:
        coords=list(part.exterior.coords)
        for a,b in zip(coords[:-1],coords[1:]):
            edge=LineString([a,b])
            if any(edge.hausdorff_distance(o)<.01 for o in openings):removed.append(edge)
            else:edges.append(edge)
        edges.extend(LineString(r.coords) for r in part.interiors)
    if len(removed)!=len(openings):raise ValueError('Reviewed opening no longer matches a complete polygon edge')
    banks=MultiLineString(edges)
    np.testing.assert_allclose(banks.length+sum(e.length for e in removed),water.boundary.length,atol=1e-5)
    return banks


def depth_field(water,banks,affine,shape_,crs,settings):
    """Fourfold sampling and metre-based smoothing, bounded-memory tiles.

    Average depth over each output pixel's water subpixels, then retain only
    output pixels whose centres belong to the polygon. Never expand onto land.
    """
    # Exact segment index avoids testing every river-bank vertex for every cell.
    lines=list(banks.geoms) if hasattr(banks,'geoms') else [banks]
    segments=np.concatenate([np.stack([np.asarray(line.coords)[:-1],np.asarray(line.coords)[1:]],axis=1) for line in lines])
    bank_index=shapely.STRtree(shapely.linestrings(segments))
    shapely.prepare(water)
    result=np.full(shape_,np.nan);factor=4;tile=128;halo=8
    spacing=np.hypot(affine.a,affine.d)/factor
    for r in range(0,shape_[0],tile):
        for c in range(0,shape_[1],tile):
            h=min(tile,shape_[0]-r);w=min(tile,shape_[1]-c)
            local=affine*Affine.translation(c,r)*Affine.scale(1/factor)*Affine.translation(-halo,-halo)
            rr,cc=np.indices((h*factor+2*halo,w*factor+2*halo))
            x,y=local*(cc+.5,rr+.5);points=shapely.points(x,y)
            inside=shapely.covers(water,points)
            if not inside.any():continue
            indices,values=bank_index.query_nearest(points[inside],all_matches=False,return_distance=True)
            nearest=np.empty(int(inside.sum()));nearest[indices[0]]=values
            distance=np.full(inside.shape,np.nan);distance[inside]=nearest
            z=profile(distance,settings['shore_transition_m'],settings['max_depth_m'],
                      settings['shore_gradient'],settings['interior_gradient'])
            s=settings['smoothing']
            z,_=smooth_creases(z,distance,spacing,s['scale_m'],s['shore_lock_m'],s['shore_blend_m'])
            z=z[halo:halo+h*factor,halo:halo+w*factor].reshape(h,factor,w,factor)
            count=np.isfinite(z).sum(axis=(1,3));total=np.nansum(z,axis=(1,3))
            result[r:r+h,c:c+w]=np.divide(total,count,out=np.full((h,w),np.nan),where=count>0)
    return result


def build(baseline_root,slug='jamestown'):
    root=terrain.BASE;label=f'{slug}-polygon-rivers'
    baseline=baseline_root/'output'/slug
    manifest_path=terrain.output_file(baseline,'manifest.json')
    manifest=json.loads(manifest_path.read_text())
    location=json.loads((root/'locations.json').read_text())[slug]
    if manifest['location']!=location:raise ValueError('Baseline location/rotation does not match')
    settings=json.loads((root/'water-configs/river-preferred.json').read_text())
    water_report=json.loads(terrain.output_file(baseline,'water-report.json').read_text())
    expected=hydrology.water_settings(terrain,slug)
    if water_report['defaults']!=expected:raise ValueError('Baseline water settings changed; rebuild baseline first')
    lake_features=json.loads(terrain.output_file(baseline,'lake-polygons.geojson').read_text())['features']
    if slug=='jamestown' and not any('waiuna' in (f['properties'].get('name') or '').lower() for f in lake_features):
        raise ValueError('Baseline is missing Waiuna Lagoon')
    snapshot=root/'river-configs'/f'{slug}-polygons.geojson'
    polygon_source=json.loads(snapshot.read_text())
    if 'location' in polygon_source and polygon_source['location']!=location:
        raise ValueError('Polygon snapshot location/rotation does not match; refresh and review openings')
    features=polygon_source['features']
    opening_path=root/'river-configs'/('jamestown-hollyford-openings.geojson' if slug=='jamestown' else f'{slug}-openings.geojson')
    opening_features=json.loads(opening_path.read_text())['features'] if opening_path.exists() else []
    crs=terrain.projection(location['longitude'],location['latitude'])
    project=Transformer.from_crs(4326,crs,always_xy=True).transform
    openings=[transform(project,shape(f['geometry'])) for f in opening_features]
    lake_open=next((g for g,f in zip(openings,opening_features) if slug=='jamestown' and 'Lake' in f['properties']['name']),None)
    sea_open=next((g for g,f in zip(openings,opening_features) if slug=='jamestown' and 'Sea' in f['properties']['name']),None)
    joins=[(lake_open,1.),(sea_open,0.)] if slug=='jamestown' else [(g,float(f['properties']['water_level_m'])) for g,f in zip(openings,opening_features)]
    destination=root/'output'/label
    staging=root/'output'/f'.{label}-building'
    if staging.exists():shutil.rmtree(staging)
    staging.mkdir(parents=True)
    arrays={};changes={};reports=[]
    for kind,width in terrain.WIDTHS.items():
        affine=terrain.grid(location,width);spacing=width/4096
        original=np.array(Image.open(terrain.output_file(baseline,f'{slug}-{kind}.png')))
        result=original.copy();owned=np.zeros(original.shape,bool)
        with rasterio.open(baseline_root/manifest['source_surface_cache']/f'surface-{kind}.tif') as src:
            raw=src.read(1)
            if src.shape!=(4096,4096) or not src.transform.almost_equals(affine):raise ValueError('Source grid mismatch')
        lakes=hydrology.feature_mask(terrain,location,kind,lake_features)
        adjoining=lakes | (raw<0)
        for feature in features:
            props=feature['properties'];name=props.get('name') or 'Unnamed river'
            water=transform(project,shape(feature['geometry']))
            if not water.is_valid:raise ValueError(f'Invalid polygon: {name}')
            window=hydrology.lake_window(water,affine,4096)
            if window is None or not window[2].any():continue
            r,c,mask=window;h,w=mask.shape
            local=affine*Affine.translation(c,r)
            is_hollyford=slug=='jamestown' and 'hollyford' in name.lower()
            selected_openings=openings if is_hollyford else [g for g,f in zip(openings,opening_features) if f['properties'].get('objectid')==props.get('OBJECTID')]
            banks=retained_banks(water,selected_openings) if selected_openings else water.boundary
            print(f'{kind}: {name}, {int(mask.sum())} water pixels',flush=True)
            depth=depth_field(water,banks,local,mask.shape,crs,settings)
            rr,cc=np.indices(mask.shape);x,y=local*(cc+.5,rr+.5)
            if is_hollyford:
                points=shapely.points(x,y)
                dl=shapely.distance(points,lake_open);ds=shapely.distance(points,sea_open)
                reference=ds/(dl+ds)  # Lake McKerrow 1 m; sea 0 m.
                reference_method='Distance-weighted 1 m lake to 0 m sea reference; synthetic, not surveyed stage'
            else:
                values=np.maximum(raw[r:r+h,c:c+w],0)
                weights=gaussian_filter(mask.astype(float),30/spacing)
                reference=gaussian_filter(np.where(mask,values,0.),30/spacing)/np.maximum(weights,1e-12)
                reference_method='30 m masked smoothing of uncarved DEM inside mapped polygon, clamped at sea level'
            target=terrain.encode(np.where(np.isfinite(depth),reference-depth,0))
            valid=mask & np.isfinite(depth)
            view=result[r:r+h,c:c+w]
            view[valid]=np.minimum(view[valid],target[valid]);owned[r:r+h,c:c+w]|=valid
            reports.append(dict(grid=kind,name=name,objectid=props.get('OBJECTID'),pixels=int(valid.sum()),
                                max_depth_m=float(depth[valid].max()),reference_min_m=float(reference[valid].min()),
                                reference_max_m=float(reference[valid].max()),reference_method=reference_method,
                                excluded_bank_edges=len(selected_openings)))
        # Join Hollyford to existing water only: a 100 m outward feather of its
        # opening cross-section, taking the lower absolute target, never adding depths.
        for opening,level in joins:
            window=hydrology.lake_window(opening.buffer(100),affine,4096)
            if window is None:continue
            r,c,mask=window;h,w=mask.shape;rr,cc=np.indices(mask.shape)
            x,y=(affine*Affine.translation(c,r))*(cc+.5,rr+.5);points=shapely.points(x,y)
            d=shapely.distance(points,opening)
            enddist=np.minimum(shapely.distance(points,Point(opening.coords[0])),shapely.distance(points,Point(opening.coords[-1])))
            # Tangential distance to the closest bank endpoint, without treating
            # the opening itself as a shore. This fades naturally at its corners.
            across=np.sqrt(np.maximum(0,enddist**2-d**2))
            depth=profile(across,settings['shore_transition_m'],settings['max_depth_m'],settings['shore_gradient'],settings['interior_gradient'])
            u=np.clip(1-d/100,0,1);blend=u*u*(3-2*u)
            view=result[r:r+h,c:c+w]
            old=view.astype(float)*terrain.STEP+terrain.HEIGHT_MIN
            target=terrain.encode(old+blend*(np.minimum(old,level-depth)-old))
            # A 1.5-pixel collar closes raster classification gaps exactly at
            # the reviewed opening (some source DEM mouth pixels are positive).
            connection=adjoining[r:r+h,c:c+w] | ((d<=1.5*spacing) & (across>spacing))
            valid=mask & connection & ~owned[r:r+h,c:c+w]
            view[valid]=np.minimum(view[valid],target[valid]);owned[r:r+h,c:c+w]|=valid
        assert np.array_equal(result[~owned],original[~owned])
        assert np.all(result<=original)
        arrays[kind]=result;changes[kind]=result!=original
        print(f'{kind}: {int(changes[kind].sum())} changed pixels',flush=True)
    # Same footprint/rotation and 4:1 nesting as the original map package.
    affected=changes['playable'].reshape(1024,4,1024,4).any(axis=(1,3))
    means=arrays['playable'].astype(float).reshape(1024,4,1024,4).mean(axis=(1,3))
    centre=arrays['world'][1536:2560,1536:2560]
    centre[affected]=np.rint(means[affected]).astype('uint16')
    np.testing.assert_array_less(np.abs(centre[affected]-means[affected]),.50001)
    for kind,array in arrays.items():
        original=np.array(Image.open(terrain.output_file(baseline,f'{slug}-{kind}.png')))
        changes[kind]=array!=original
        path=staging/f'{label}-{kind}.png';Image.fromarray(array).save(path)
        check=np.array(Image.open(path));assert check.dtype==np.uint16 and check.shape==(4096,4096)
        np.testing.assert_array_equal(check,array)
        terrain.preview(array,location,kind,staging/f'{kind}-preview.jpg',
                        lake_mask=hydrology.feature_mask(terrain,location,kind,lake_features))
    for filename in ['water-report.json','lake-polygons.geojson','manifest.json']:
        shutil.copy(terrain.output_file(baseline,filename),staging/('baseline-manifest.json' if filename=='manifest.json' else filename))
    for path in [snapshot,opening_path,root/'water-configs/river-preferred.json',root/'river-configs/SOURCES.md']:
        if path.exists():shutil.copy(path,staging/path.name)
    report=dict(generated_utc=datetime.now(timezone.utc).isoformat(),location=location,settings=settings,
                baseline_manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
                baseline_heightmap_sha256={k:hashlib.sha256(terrain.output_file(baseline,f'{slug}-{k}.png').read_bytes()).hexdigest() for k in terrain.WIDTHS},
                polygon_sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest(),
                opening_count=len(opening_features),opening_sha256=hashlib.sha256(opening_path.read_bytes()).hexdigest() if opening_path.exists() else None,rivers=reports,
                changed_pixels={k:int(v.sum()) for k,v in changes.items()},composition='Lower absolute target wins; no additive excavation',
                limitations=['Synthetic beds; no hydraulic/in-game validation or automatic water sources.',
                             'Mapped river polygons only; no centreline extensions. Snapshot covers the playable-area query, not every river in the surrounding world.',
                             'Only explicitly recorded opening edges are excluded. Other polygon cutoffs remain pending review.',
                             'Other river reference levels are DEM estimates; no downstream monotonicity or cliff filter is enforced.',
                             'Existing depression detection is deferred. Wide reaches can have flat beds.',
                             'Recorded connections feather into existing lake/sea beds over 100 m inside existing water masks, plus a 1.5-pixel collar at each reviewed opening to prevent raster seams.'])
    (staging/'river-report.json').write_text(json.dumps(report,indent=2)+'\n')
    instructions=(baseline/'IMPORT.txt').read_text().replace(f'{slug}-playable.png',f'{label}-playable.png').replace(f'{slug}-world.png',f'{label}-world.png')
    instructions+='\nPOLYGON RIVER TRIAL\nUse BOTH PNGs from this package. River depth cap 10 m; 2 m shoreline transition; 1 m crease smoothing.\nLakes, lagoons and seabed inherit the existing synthetic baseline. Recorded opening edges have local connections. See the opening snapshot and river report.\nPlace water sources manually; depth parameters are synthetic bed depths, not guaranteed simulated water depths.\nRiver geometry: LINZ River Name Polygons Pilot, CC BY 4.0. See supplied source metadata and river-report.json.\n'
    (staging/'IMPORT.txt').write_text(instructions)
    archive=terrain.package_outputs(staging,label)
    if destination.exists():shutil.rmtree(destination)
    staging.replace(destination)
    print(f'COMPLETE: {destination/archive.name}',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-root',type=Path,required=True)
    parser.add_argument('--location',default='jamestown')
    args=parser.parse_args();build(args.baseline_root.resolve(),args.location)


if __name__=='__main__':main()
