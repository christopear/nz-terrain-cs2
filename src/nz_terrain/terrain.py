"""Build true-scale vanilla CS:II terrain packages from named NZ locations."""
from __future__ import annotations
import argparse
import sys
from .paths import project_root, preview_font_path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import time
from tempfile import NamedTemporaryFile
from urllib.parse import urljoin
from urllib.request import urlopen
import zipfile

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from pyproj import CRS, Transformer
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import Affine, from_bounds
from rasterio.vrt import WarpedVRT
from scipy.ndimage import distance_transform_edt, gaussian_filter, map_coordinates

BASE = project_root()
ROOT = 'https://nz-elevation.s3.ap-southeast-2.amazonaws.com/'
CONTOUR = 'new-zealand/new-zealand-contour/dem_8m/2193/'
WIDTHS = {'playable': 14336, 'world': 57344}
PIXELS = 4096
SEA = 8192
HEIGHT_MIN, HEIGHT_SPAN = -512, 4096
STEP = HEIGHT_SPAN / 65535
VERSION = 'locations-v1'
CACHE = BASE / 'data' / 'catalogue'
LETTERS = 'ABCDEFGHJKLMNPQRSTUVWXYZ'


def projection(lon, lat):
    return CRS.from_proj4(f'+proj=tmerc +lat_0={lat} +lon_0={lon} +k=1 +datum=WGS84 +units=m')


def grid(config, width, size=PIXELS):
    # Map-up points clockwise from true north. Rotate the sampling footprint,
    # never the finished image (which would crop/resize the real-world terrain).
    return Affine.rotation(-config['rotation_deg']) * from_bounds(-width/2, -width/2, width/2, width/2, size, size)


def outline(config, width, target=4326):
    crs = projection(config['longitude'], config['latitude'])
    t = Transformer.from_crs(crs, target, always_xy=True)
    a = grid(config, width)
    corners = [(0, 0), (PIXELS, 0), (PIXELS, PIXELS), (0, PIXELS), (0, 0)]
    points = []
    for first, last in zip(corners, corners[1:]):
        for f in np.linspace(0, 1, 33, endpoint=False):
            x, y = a * (first[0]+f*(last[0]-first[0]), first[1]+f*(last[1]-first[1]))
            points.append(t.transform(x, y))
    return points + [points[0]]


def bounds(points):
    a = np.array(points)
    return [float(a[:, 0].min()), float(a[:, 1].min()), float(a[:, 0].max()), float(a[:, 1].max())]


def intersects(a, b):
    return a[0] < b[2] and a[2] > b[0] and a[1] < b[3] and a[3] > b[1]


def fetch(url):
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / (hashlib.sha256(url.encode()).hexdigest()+'.json')
    if dest.exists():
        return json.loads(dest.read_text())
    # Reuse metadata collected for the initial Jamestown build.
    legacy = BASE / 'data' / 'metadata' / url.replace('https://', '').replace('/', '_')
    if len(legacy.name.encode()) <= 255 and legacy.exists():
        data = legacy.read_bytes()
    else:
        for attempt in range(4):
            try:
                with urlopen(url, timeout=90) as response:
                    data = response.read()
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(2**attempt)
    result = json.loads(data)
    temp = dest.with_suffix('.tmp')
    temp.write_bytes(data)
    temp.replace(dest)
    return result


def sheet_bounds(filename):
    """LINZ Topo50 sheet grid; use only as coarse metadata prefilter.

    CA10 is [1228000,5046000,1252000,5082000] in EPSG:2193.
    Precise STAC geometries are checked after this inexpensive name filter.
    Unknown file conventions are fetched rather than silently discarded.
    """
    match = re.fullmatch(r'([A-Z]{2})(\d{2})(?:_10000_(\d{2})(\d{2}))?\.json', filename)
    if not match:
        return None
    row, col, subrow, subcol = match.groups()
    if any(c not in LETTERS for c in row):
        return None
    index = lambda s: LETTERS.index(s[0])*len(LETTERS)+LETTERS.index(s[1])
    west = 1228000+(int(col)-10)*24000
    north = 5082000-(index(row)-index('CA'))*36000
    if subrow:
        west += (int(subcol)-1)*4800
        north -= (int(subrow)-1)*7200
        return [west, north-7200, west+4800, north]
    return [west, north-36000, west+24000, north]


def discover_collection(folder, bbox, nztm_bbox, workers):
    base = ROOT+folder
    collection = fetch(base+'collection.json')
    if not any(intersects(b, bbox) for b in collection['extent']['spatial']['bbox']):
        return []
    selected = []
    for link in collection['links']:
        if link['rel'] != 'item':
            continue
        coarse = sheet_bounds(link['href'].split('/')[-1])
        if coarse is None or intersects(coarse, nztm_bbox):
            selected.append(urljoin(base, link['href']))
    def item(url):
        metadata = fetch(url)
        if not intersects(metadata['bbox'], bbox):
            return None
        asset = next(v for v in metadata['assets'].values() if 'geotiff' in v.get('type', ''))
        return {'url': urljoin(url, asset['href']), 'bbox': metadata['bbox'],
                'collection': collection.get('title', collection['id']),
                'item_url': url, 'source_checksum': asset.get('file:checksum')}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return [v for v in pool.map(item, selected) if v is not None]


def discover(config, workers):
    bbox = bounds(outline(config, WIDTHS['world']))
    nztm_bbox = bounds(outline(config, WIDTHS['world'], 2193))
    result = {'contour': discover_collection(CONTOUR, bbox, nztm_bbox, workers), 'lidar': [], 'copernicus': []}
    # Later collections take precedence, explicitly configured oldest to newest.
    for folder in config.get('lidar_collections', []):
        result['lidar'].extend(discover_collection(folder, bbox, nztm_bbox, workers))
    for lat in range(math.floor(bbox[1]), math.floor(bbox[3])+1):
        for lon in range(math.floor(bbox[0]), math.floor(bbox[2])+1):
            ns, ew = ('N' if lat >= 0 else 'S'), ('E' if lon >= 0 else 'W')
            stem = f'Copernicus_DSM_COG_10_{ns}{abs(lat):02d}_00_{ew}{abs(lon):03d}_00_DEM'
            result['copernicus'].append({'url': f'https://copernicus-dem-30m.s3.amazonaws.com/{stem}/{stem}.tif',
                                         'bbox': [lon, lat, lon+1, lat+1], 'collection': 'Copernicus GLO-30 Public'})
    return result


def write_float(path, z, config, width):
    temp = path.with_suffix('.partial.tif')
    with rasterio.open(temp, 'w', driver='GTiff', width=PIXELS, height=PIXELS, count=1,
                       dtype='float32', crs=projection(config['longitude'], config['latitude']),
                       transform=grid(config, width), nodata=float('nan'), compress='deflate', tiled=True) as dst:
        dst.write(z.astype('float32'), 1)
    temp.replace(path)


def verified_ocean_tile(item):
    """Only treat absent tiles as ocean if absent from BOTH public catalogues.

    GLO-30 alone also omits some restricted land tiles; GLO-90 is worldwide.
    A listed tile that returns 404 remains an error, never an ocean substitute.
    """
    stem=item['url'].split('/')[-2]
    sources=[]
    for resolution,code in [(30,'10'),(90,'30')]:
        url=f'https://copernicus-dem-{resolution}m.s3.amazonaws.com/tileList.txt'
        path=CACHE/f'copernicus-{resolution}m-tileList.txt'
        if not path.exists():
            payload=urlopen(url,timeout=90).read()
            if not payload.startswith(b'Copernicus_DSM_COG_'):raise ValueError('Invalid Copernicus tile catalogue')
            with NamedTemporaryFile(dir=CACHE,delete=False) as temporary:
                temporary.write(payload)
            Path(temporary.name).replace(path)
        names=set(path.read_text().splitlines())
        if stem.replace('_COG_10_','_COG_'+code+'_') in names:return False
        sources.append(url)
    item['ocean_tile_verified_against']=sources
    return True


def ocean_tile_grid(item,config,width,size=PIXELS):
    result=np.full((size,size),np.nan,dtype='float32')
    affine=grid(config,width,size)
    inverse=Transformer.from_crs(projection(config['longitude'],config['latitude']),4326,always_xy=True)
    west,south,east,north=item['bbox']
    for start in range(0,size,256):
        yy,xx=np.mgrid[start:min(start+256,size),:size]
        lon,lat=inverse.transform(*(affine*(xx+.5,yy+.5)))
        inside=(lon>=west)&(lon<=east)&(lat>=south)&(lat<=north)
        result[start:start+256][inside]=0
    return result


def prepare(config, directory, sources, workers):
    os.environ['GDAL_DISABLE_READDIR_ON_OPEN'] = 'EMPTY_DIR'
    os.environ['GDAL_HTTP_TIMEOUT'] = '120'
    os.environ['GDAL_HTTP_MAX_RETRY'] = '4'
    crs = projection(config['longitude'], config['latitude'])
    for name, width in WIDTHS.items():
        bbox = bounds(outline(config, width))
        for kind, items in sources.items():
            path = directory / f'{kind}-{name}.tif'
            if path.exists():
                continue
            relevant = [i for i in items if intersects(i['bbox'], bbox)]
            z = np.full((PIXELS, PIXELS), np.nan, dtype='float32')
            def read_tile(item):
                url = item['url']
                try:
                    header=rasterio.open(url)
                except rasterio.errors.RasterioIOError as error:
                    if kind=='copernicus' and '404' in str(error) and verified_ocean_tile(item):
                        print(f'Catalogue-confirmed ocean tile: {item["bbox"]}',flush=True)
                        return ocean_tile_grid(item,config,width)
                    raise
                with header:
                    factors = header.overviews(1)
                    # Published LINZ 1m grids: use 2m overviews for 3.5m output,
                    # 8m for 14m output. Other source grids retain native detail.
                    level = -1
                    if kind == 'lidar':
                        candidates = [i for i, f in enumerate(factors) if f <= (2 if name == 'playable' else 8)]
                        level = max(candidates, default=-1)
                options = {'overview_level': level} if level >= 0 else {}
                with rasterio.open(url, **options) as src:
                    with WarpedVRT(src, crs=crs, transform=grid(config, width), width=PIXELS, height=PIXELS,
                                   dtype='float32', nodata=float('nan'), resampling=Resampling.bilinear) as vrt:
                        return vrt.read(1, masked=True).filled(np.nan)
            # Bounded batches retain deterministic source priority without
            # keeping an entire region's arrays resident in memory.
            with ThreadPoolExecutor(max_workers=workers) as pool:
                for start in range(0, len(relevant), workers):
                    for a in pool.map(read_tile, relevant[start:start+workers]):
                        valid = np.isfinite(a)
                        z[valid] = a[valid]
                    print(f"{config['name']} {name}/{kind}: {min(start+workers,len(relevant))}/{len(relevant)} tiles", flush=True)
            write_float(path, z, config, width)


def encode(z):
    if not np.isfinite(z).all() or z.min() < HEIGHT_MIN or z.max() > HEIGHT_MIN+HEIGHT_SPAN:
        raise ValueError('Missing elevations or elevation outside encoding range; clipping refused')
    return np.rint((z-HEIGHT_MIN)/HEIGHT_SPAN*65535).astype('uint16')


def rebuild_seabed(raw, lake_masks=None, maximum=30., gradient=.05):
    """Final operation: never blend survey data over synthetic seabed again."""
    masks = {n: a < SEA for n, a in raw.items()}
    if lake_masks is not None:
        masks = {n: mask & ~lake_masks[n] for n,mask in masks.items()}
    result = {n: a.copy() for n, a in raw.items()}
    pad = 320
    coords = (np.arange(PIXELS+pad*2)-pad+.5)/4+1536-.5
    yy, xx = np.meshgrid(coords, coords, indexing='ij')
    expanded = map_coordinates(masks['world'].astype('uint8'), [yy, xx], order=0, mode='nearest').astype(bool)
    expanded[pad:pad+PIXELS, pad:pad+PIXELS] = masks['playable']
    del yy, xx
    for name, spacing in [('playable', 3.5), ('world', 14)]:
        mask = expanded if name == 'playable' else masks[name]
        if not mask.any():
            continue
        if mask.all():
            depth = np.full(mask.shape, maximum, dtype='float32')
        else:
            distance = distance_transform_edt(mask, sampling=spacing).astype('float32')
            distance = gaussian_filter(distance, sigma=14/spacing, mode='nearest')
            depth = maximum*np.tanh(distance/(maximum/gradient))
        if name == 'playable':
            depth = depth[pad:pad+PIXELS, pad:pad+PIXELS]
        values = np.minimum(encode(-depth), SEA-1)
        result[name][masks[name]] = values[masks[name]]
    centre = result['world'][1536:2560, 1536:2560]
    water = masks['world'][1536:2560, 1536:2560]
    means = result['playable'].astype('float32').reshape(1024, 4, 1024, 4).mean(axis=(1, 3))
    centre[water] = np.minimum(np.rint(means[water]), SEA-1).astype('uint16')
    for name in result:
        assert np.array_equal(result[name][~masks[name]], raw[name][~masks[name]])
        assert np.array_equal(result[name]<SEA, raw[name]<SEA)
    pure = masks['playable'].reshape(1024, 4, 1024, 4).all(axis=(1, 3)) & water
    err = float(np.max(np.abs(centre[pure].astype('float32')-means[pure]))) if pure.any() else 0
    assert err <= .5
    return result, err*STEP


def preview(a, config, name, path, lake_mask=None):
    z = a.astype('float32')*STEP+HEIGHT_MIN
    small = z.reshape(1024, 4, 1024, 4).mean(axis=(1, 3))
    dy, dx = np.gradient(small, WIDTHS[name]/1024)
    light = np.clip((1-.7*dx+.7*dy)/np.sqrt(1+dx*dx+dy*dy), 0, 1)
    stops = [-30, 0, 100, 500, 1000, 1800, 2800]
    colours = np.array([[30,75,115],[77,117,92],[108,139,94],[147,153,109],[166,156,133],[193,190,181],[240,240,235]])
    rgb = np.stack([np.interp(small, stops, colours[:,i]) for i in range(3)], axis=-1)
    rgb *= .5+.5*light[:,:,None]
    water = small < 0
    if lake_mask is None and (path.parent/'lake-polygons.geojson').exists():
        from . import hydrology
        lake_mask=hydrology.feature_mask(sys.modules[__name__],config,name,
            json.loads((path.parent/'lake-polygons.geojson').read_text())['features'])
    if lake_mask is not None:
        water |= lake_mask.reshape(1024,4,1024,4).any(axis=(1,3))
    f = np.clip(-small/30, 0, 1)
    rgb[water] = np.stack([40-15*f, 145-70*f, 180-65*f], axis=-1)[water]
    im = Image.fromarray(rgb.astype('uint8'))
    draw = ImageDraw.Draw(im)
    font = ImageFont.truetype(preview_font_path(), 18)
    if name == 'world':
        draw.rectangle((384,384,640,640), outline='white', width=2)
    draw.line((502,512,522,512), fill='white'); draw.line((512,502,512,522), fill='white')
    draw.rectangle((0,0,1024,45),fill=(15,28,35))
    draw.text((14,12),f"{config['name'].upper()} | {name} | {WIDTHS[name]/1000:g} km | map-up {config['rotation_deg']:g} deg CW from N",font=font,fill='white')
    draw.rectangle((0,984,1024,1024),fill=(15,28,35))
    draw.text((14,994),f"{config['longitude']:.5f}, {config['latitude']:.5f} | synthetic water beds",font=font,fill='white')
    im.save(path, quality=93)


def footprints(config, out):
    features = [{'type':'Feature','properties':{'name':name,'width_m':w,'rotation_deg':config['rotation_deg']},
                 'geometry':{'type':'Polygon','coordinates':[outline(config,w)]}} for name,w in WIDTHS.items()]
    (out/'footprints.geojson').write_text(json.dumps({'type':'FeatureCollection','features':features},indent=2))


def finish(config, directory, sources, out, slug, water_options=None):
    from . import hydrology
    arrays, stats = {}, {}
    for name in WIDTHS:
        data = {}
        for kind in sources:
            with rasterio.open(directory/f'{kind}-{name}.tif') as src:
                data[kind] = src.read(1)
        c, l, satellite = data['contour'], data['lidar'], data['copernicus']
        if not np.isfinite(satellite).all():
            raise ValueError('Incomplete Copernicus coverage: no silent missing-tile ocean fill')
        missing = ~np.isfinite(c)
        base = np.where(missing, satellite, c)
        # Seed offshore areas from the prior contour coast where available;
        # where contour coverage is absent, use Copernicus zero-elevation water.
        ocean = missing & (satellite <= 1)
        base[ocean] = -1
        valid = np.isfinite(l)
        if valid.any():
            weight = np.minimum(distance_transform_edt(valid, sampling=WIDTHS[name]/PIXELS)/250, 1).astype('float32')
            base[valid] = base[valid]*(1-weight[valid])+l[valid]*weight[valid]
        # Keep the inferred ocean mask, but preserve LiDAR dry-land corrections.
        base[ocean & ~valid] = -1
        arrays[name] = base
        stats[name] = {'lidar_coverage_percent':float(valid.mean()*100),
                       'contour_coverage_percent':float((~missing).mean()*100),
                       'satellite_land_fallback_percent':float((missing & ~ocean & ~valid).mean()*100)}
    arrays['world'][1536:2560,1536:2560] = arrays['playable'].reshape(1024,4,1024,4).mean(axis=(1,3))
    # Keep the source surface immutable for independently calculated river beds.
    for name,z in arrays.items():
        write_float(directory/f'surface-{name}.tif',z,config,WIDTHS[name])
    lake_targets,lake_masks,water_report,lake_features=hydrology.prepare_lakes(
        sys.modules[__name__],slug,config,arrays,water_options)
    raw = {n:encode(z) for n,z in arrays.items()}
    del arrays
    settings=water_report['defaults'];sea_depth=settings['sea_depth_m'];sea_gradient=settings['sea_gradient']
    repaired, error = rebuild_seabed(raw,lake_masks,sea_depth,sea_gradient)
    for name in WIDTHS:
        repaired[name]=hydrology.combine_targets(repaired[name],[(encode(lake_targets[name]),lake_masks[name])])
    affected=lake_masks['playable'].reshape(1024,4,1024,4).any(axis=(1,3))
    means=repaired['playable'].astype('float64').reshape(1024,4,1024,4).mean(axis=(1,3))
    centre=repaired['world'][1536:2560,1536:2560]
    centre[affected]=np.rint(means[affected]).astype('uint16')
    if affected.any():assert np.max(np.abs(centre[affected]-means[affected]))<=.5
    (out/'water-report.json').write_text(json.dumps(water_report,indent=2))
    (out/'lake-polygons.geojson').write_text(json.dumps(lake_features,indent=2))
    hydrology.plot_water_profiles(sys.modules[__name__],config,raw,repaired,water_report,lake_features,out)
    report = {'pipeline_version':VERSION,'generated_utc':datetime.now(timezone.utc).isoformat(),
              'location':config,'widths_m':WIDTHS,'pixels':PIXELS,'sources':sources,
              'water_pipeline_version':'independent-beds-v1','lake_count':len(water_report['lakes']),
              'source_surface_cache':str(directory.relative_to(BASE)),
              'black_elevation_m':HEIGHT_MIN,'white_elevation_m':HEIGHT_MIN+HEIGHT_SPAN,
              'height_span_m':HEIGHT_SPAN,'sea_level_game_m_if_zero_based_import':512,
              'seabed':{'formula':f'depth = {sea_depth} * tanh(smoothed distance to shore / {sea_depth/sea_gradient})',
                        'distance_smoothing_sigma_m':14,'maximum_depth_m':sea_depth,'near_shore_gradient':sea_gradient,
                        'land_pixels_changed':0,'shoreline_mask_pixels_changed':0,
                        'fully_submerged_centre_alignment_error_m':error},'grids':stats,
              'limitations':['Synthetic seabed, not bathymetric survey.',
                             'No formal vertical datum harmonisation across source datasets.',
                             'Lake beds are synthetic; water levels and sources need editor testing. Rivers require the separate configured river build.',
                             'No in-game calibration test. Small projection distortion increases away from centre.']}
    for name,a in repaired.items():
        target = out/f'{slug}-{name}.png'
        Image.fromarray(a).save(target)
        with Image.open(target) as im:
            check = np.array(im)
            assert check.dtype == np.uint16 and check.shape == (PIXELS,PIXELS) and np.array_equal(a,check)
        stats[name].update({'filename':target.name,'min_m':float(a.min()*STEP+HEIGHT_MIN),'max_m':float(a.max()*STEP+HEIGHT_MIN),
                           'sha256':hashlib.sha256(target.read_bytes()).hexdigest()})
        preview(a,config,name,out/f'{name}-preview.jpg',lake_masks[name])
    (out/'manifest.json').write_text(json.dumps(report,indent=2))
    footprints(config,out)
    instructions = f'''{config['name'].upper()} — VANILLA CITIES: SKYLINES II HEIGHTMAPS
Centre (longitude, latitude): {config['longitude']}, {config['latitude']}
Map-up direction: {config['rotation_deg']} degrees CLOCKWISE from true north.
Rotation applies to the real-world sampling footprint, with no resizing/cropping.

Import {slug}-playable.png (14.336 km square) and {slug}-world.png (57.344 km square) together.
Both are 4096 x 4096, unsigned 16-bit grayscale. Do not auto-contrast or resize.
Place them in %USERPROFILE%\\AppData\\LocalLow\\Colossal Order\\Cities Skylines II\\Heightmaps.
The preview images and GeoJSON are not import heightmaps.

HEIGHT SETTINGS
Encoding: real metres = pixel / 65535 * 4096 - 512.
Use the same settings as the successful Jamestown import: height span 4096 m;
if black imports as game height 0, sea level is game height 512 m.
Verify your editor's height scale/offset before editing. No in-game test here.
These packages contain terrain, not complete playable cities/maps.

SEABED
Synthetic sea depth ceiling: {sea_depth:g} m; near-shore gradient parameter: {sea_gradient:g} m/m.
Below-sea pixels outside mapped lakes are rebuilt AFTER source blending to avoid
submerged tile steps. Land is unchanged by the seabed step.
The world centre is aligned with the playable terrain. Mixed coastal cells
prioritise preserving land over exact averaging, except cells updated for lake beds.
Depths are not real surveyed bathymetry.

LAKES
Mapped perennial lakes have synthetic beds calculated from a separate water level.
Default lake depth parameter: {settings['lake_depth_m']:g} m; lagoon depth: {settings['lagoon_depth_m']:g} m.
Sourced per-lake overrides and effective gradients are in water-report.json.
Small lakes may be shallower because depth grows gradually away from the shore.
Lake interiors are excluded from the ocean calculation, even below sea level.
Lake targets replace the source surface; depths are never subtracted repeatedly.
See metadata/water-report.json for lake water levels, maximum depths and caveats.
Place water sources in the editor and check outlet levels; terrain alone does not
fill lakes. River candidates use the original source surface, then take the lower
absolute bed in overlaps. Unmapped water and other hydrographic classes need review.

DATA
Prefer bare-earth LINZ LiDAR where configured/available; feather outer coverage
over 250 m into LINZ 8 m contour terrain, using Copernicus where necessary.
Copernicus includes vegetation. Cross-source vertical datums are not formally
harmonised; this is a cartographic/game terrain, not a survey product.
See manifest.json for actual coverage, source URLs, configuration and checksums.

ATTRIBUTION
Contains data sourced from Toitu Te Whenua Land Information New Zealand,
licensed for reuse under CC BY 4.0: https://creativecommons.org/licenses/by/4.0/
NZ 8 m DEM created by Geographx from LINZ Topo50 contours. LiDAR is provided by
its respective councils/survey providers via LINZ; see linked collection metadata.
https://github.com/linz/elevation
Produced using Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and
© Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS
by the European Union and ESA; all rights reserved.
https://registry.opendata.aws/copernicus-dem/
'''
    if (BASE/'water-configs'/f'{slug}-polygons.geojson').exists():
        instructions+='\nSupplementary water polygons: © OpenStreetMap contributors, ODbL 1.0.\nhttps://www.openstreetmap.org/copyright\nThe adapted polygon database is supplied in metadata/lake-polygons.geojson.\n'
    (out/'IMPORT.txt').write_text(instructions)


def build(slug, config, workers, command, water_options=None):
    digest = hashlib.sha256(json.dumps([VERSION,config],sort_keys=True).encode()).hexdigest()[:12]
    directory = BASE/'data'/'locations'/slug/digest
    directory.mkdir(parents=True,exist_ok=True)
    destination = BASE/'output'/slug
    if command == 'footprints':
        metadata = destination/'metadata'
        metadata.mkdir(parents=True,exist_ok=True)
        footprints(config,metadata)
        return
    sourcefile = directory/'sources.json'
    if sourcefile.exists():
        sources = json.loads(sourcefile.read_text())
    else:
        sources = discover(config,max(workers,8))
        sourcefile.write_text(json.dumps(sources,indent=2))
    print(slug, {k:len(v) for k,v in sources.items()},flush=True)
    prepare(config,directory,sources,workers)
    sourcefile.write_text(json.dumps(sources,indent=2))
    staging = BASE/'output'/f'.{slug}-building'
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True,exist_ok=True)
    finish(config,directory,sources,staging,slug,water_options)
    archive = package_outputs(staging,slug)
    # Only replace this location after successful generation and verification.
    if destination.exists():
        shutil.rmtree(destination)
    staging.rename(destination)
    print(f'COMPLETE: {destination/ archive.name}',flush=True)


def output_file(directory, filename):
    """Resolve a file in a structured package (or a legacy flat package)."""
    folder = output_folder(filename)
    structured = directory/folder/filename
    return structured if structured.exists() else directory/filename


def output_folder(filename):
    if filename.endswith(('-playable.png','-world.png')):
        return 'heightmaps'
    if Path(filename).suffix.lower() in {'.jpg','.png'}:
        return 'previews'
    return 'metadata'


def package_outputs(directory, slug):
    """Organise generated files and publish a ZIP with the same folder layout."""
    for source in sorted(directory.iterdir()):
        if not source.is_file() or source.suffix=='.zip' or source.name=='IMPORT.txt':
            continue
        target=directory/output_folder(source.name)/source.name
        target.parent.mkdir(exist_ok=True)
        if target.exists():
            raise FileExistsError(f'Refusing to overwrite duplicate package file: {target}')
        source.rename(target)
    manifest=directory/'metadata'/'manifest.json'
    if manifest.exists():
        data=json.loads(manifest.read_text())
        for grid_info in data.get('grids',{}).values():
            if 'filename' in grid_info:
                grid_info['filename']='heightmaps/'+Path(grid_info['filename']).name
        manifest.write_text(json.dumps(data,indent=2))
    instructions=directory/'IMPORT.txt'
    layout='PACKAGE FOLDERS\nHeightmaps to import are in heightmaps/.\nPreview images are in previews/; reports, profiles and geographic data are in metadata/.\nCopy only the two heightmap PNGs into the game\'s Heightmaps folder.\n\n'
    if instructions.exists() and not instructions.read_text().startswith('PACKAGE FOLDERS'):
        instructions.write_text(layout+instructions.read_text())
    archive=directory/f'{slug}-heightmaps.zip'
    temporary=directory/f'.{slug}-heightmaps.zip.tmp'
    files=sorted(p for p in directory.rglob('*') if p.is_file() and p!=temporary and p.suffix!='.zip')
    with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for p in files:
            z.write(p,arcname=p.relative_to(directory).as_posix())
    with zipfile.ZipFile(temporary) as z:
        assert z.testzip() is None
        assert set(z.namelist())=={p.relative_to(directory).as_posix() for p in files}
    temporary.replace(archive)
    return archive


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',nargs='?',choices=['build','footprints','list'],default='build')
    parser.add_argument('--location',action='append',help='Preset slug; repeat for multiple locations')
    parser.add_argument('--all',action='store_true')
    parser.add_argument('--with-rivers',action='store_true',help='Also rebuild configured river candidates from the new water-bed baseline')
    parser.add_argument('--workers',type=int,default=3,choices=range(1,9))
    for kind in ['sea','lake','lagoon']:
        parser.add_argument(f'--default-{kind}-depth',type=float,help='Synthetic depth ceiling in metres')
        parser.add_argument(f'--default-{kind}-gradient',type=float,help='Near-shore slope in m/m, e.g. 0.02 = 2%%')
    args = parser.parse_args()
    if args.with_rivers and args.command!='build':parser.error('--with-rivers requires the build command')
    water_options={}
    for kind in ['sea','lake','lagoon']:
        for parameter in ['depth','gradient']:
            value=getattr(args,f'default_{kind}_{parameter}')
            if value is not None:
                if not math.isfinite(value) or value<=0:parser.error('Water depths and gradients must be finite and positive')
                water_options[f'{kind}_{"depth_m" if parameter=="depth" else parameter}']=value
    locations = json.loads((BASE/'locations.json').read_text())
    if args.command == 'list':
        for slug,c in locations.items():
            print(f"{slug:20s} {c['longitude']:10.5f}, {c['latitude']:10.5f}  {c['rotation_deg']:g} degrees")
        return
    if args.all and args.location:
        parser.error('Use --all or --location, not both')
    selected = list(locations) if args.all else args.location
    if not selected or any(s not in locations for s in selected):
        parser.error('Choose --location from '+', '.join(locations)+' or --all')
    for slug in selected:
        config = locations[slug]
        if not (-180<=config['longitude']<=180 and -80<=config['latitude']<=80 and math.isfinite(config['rotation_deg'])):
            parser.error(f'Invalid coordinates/rotation for {slug}')
        build(slug,config,args.workers,args.command,water_options)
        if args.with_rivers:
            if (BASE/'river-configs'/f'{slug}.json').exists():
                from . import rivers
                rivers.build(slug)
            else:print(f'{slug}: no configured river routes; sea/lake/lagoon package built',flush=True)

if __name__ == '__main__':
    main()
