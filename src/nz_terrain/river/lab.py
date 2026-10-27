"""Jamestown/lower-Hollyford diagnostics; never carves or rebuilds water beds."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from pyproj import Transformer
import rasterio
from shapely.geometry import Polygon, Point, shape, mapping
from shapely.ops import transform, unary_union

from .. import terrain
from .geometry import line_parts, polygon_width, sample_surface, transect

SERVICE = 'https://services.arcgis.com/xdsHIIxuCWByZiCB/ArcGIS/rest/services/LINZ_NZ_River_Name_Polygons_Pilot/FeatureServer/0'
SNAPSHOT = terrain.BASE / 'river-configs/jamestown-hollyford-polygons.geojson'


def fetch_polygons():
    """Explicit acquisition step; regular inspection works entirely offline."""
    location = json.loads((terrain.BASE/'locations.json').read_text())['jamestown']
    params = dict(f='geojson', where='1=1', geometry=','.join(map(str, terrain.bounds(terrain.outline(location,14336)))),
                  geometryType='esriGeometryEnvelope', inSR=4326, outSR=4326,
                  outFields='*', resultRecordCount=1000, orderByFields='OBJECTID')
    features = []; seen = set(); urls = []
    while True:
        params['resultOffset'] = len(seen)
        url = SERVICE+'/query?'+urlencode(params); urls.append(url)
        with urlopen(url, timeout=90) as response:
            page = json.load(response)
        if 'features' not in page:
            raise ValueError(f'River polygon query failed: {page.get("error", page)}')
        batch = page['features']
        for feature in batch:
            key = feature['properties']['OBJECTID']
            if key in seen: raise ValueError('Polygon service repeated a page')
            seen.add(key)
            if 'hollyford' in (feature['properties'].get('name') or '').lower():
                features.append(feature)
        more = page.get('exceededTransferLimit') or page.get('properties',{}).get('exceededTransferLimit')
        if not batch and more: raise ValueError('Polygon pagination stalled')
        if not more and len(batch)<1000: break
    if not features: raise ValueError('No Hollyford polygons found; do not substitute inferred widths')
    with urlopen(SERVICE+'?f=json', timeout=90) as response:
        metadata = json.load(response)
    snapshot = dict(type='FeatureCollection', features=features, source=SERVICE, query_urls=urls,
                    retrieved_utc=datetime.now(timezone.utc).isoformat(), location=location,
                    service_geometry={k:metadata.get(k) for k in ['geometryType','hasZ','hasM','fields']},
                    attribution='LINZ River Name Polygons Pilot, CC BY 4.0. Mapped area, not surveyed bathymetry or a fixed-stage water boundary.')
    SNAPSHOT.write_text(json.dumps(snapshot,ensure_ascii=False,indent=2)+'\n')
    print(f'Saved {len(features)} Hollyford polygons to {SNAPSHOT}')


def inspect(baseline_root, output, spacing=50., half_width=500., maximum_slope=45.):
    location = json.loads((terrain.BASE/'locations.json').read_text())['jamestown']
    baseline = baseline_root/'output/jamestown'
    manifest_path = terrain.output_file(baseline, 'manifest.json')
    manifest = json.loads(manifest_path.read_text())
    if manifest['location'] != location:
        raise ValueError('Jamestown baseline centre/rotation/sources differ from the current preset')
    surface_path = (baseline_root/manifest['source_surface_cache']/'surface-playable.tif').resolve()
    with rasterio.open(surface_path) as src:
        expected_crs = terrain.projection(location['longitude'],location['latitude'])
        if src.shape != (4096,4096) or src.crs != expected_crs or not src.transform.almost_equals(terrain.grid(location,14336)):
            raise ValueError('Original surface grid does not match the Jamestown preset')
        surface = src.read(1); affine = src.transform
    project = Transformer.from_crs(4326,expected_crs,always_xy=True).transform
    unproject = Transformer.from_crs(expected_crs,4326,always_xy=True).transform
    footprint = Polygon([affine*p for p in [(0,0),(4096,0),(4096,4096),(0,4096)]])
    snapshot = json.loads(SNAPSHOT.read_text())
    if snapshot['location'] != location: raise ValueError('Refresh polygon snapshot for the current preset')
    polygon = unary_union([transform(project,shape(f['geometry'])) for f in snapshot['features']])
    if polygon.is_empty or not polygon.is_valid: raise ValueError('Invalid river polygon snapshot')
    # Use existing lake boundaries only to exclude interiors; no lake generation.
    lake_path = terrain.output_file(baseline,'lake-polygons.geojson')
    lake_features = json.loads(lake_path.read_text())['features']
    lakes = unary_union([transform(project,shape(f['geometry'])) for f in lake_features])
    routes_path = terrain.BASE/'river-configs/jamestown-network.geojson'
    routes = json.loads(routes_path.read_text())['features']
    lines = [line for f in routes if 'hollyford' in (f['properties'].get('name') or '').lower()
             for line in line_parts(transform(project,shape(f['geometry'])).intersection(footprint).difference(lakes)) if line.length>=spacing]
    if not lines: raise ValueError('No lower-Hollyford route segments in the playable footprint')
    stations=[]; sections=[]; cross_rows=[]
    offsets=np.linspace(-half_width,half_width,int(np.ceil(2*half_width/5))+1)
    for segment,line in enumerate(lines,1):
        distances=np.linspace(0,line.length,max(3,int(np.ceil(line.length/spacing))+1))
        xy=np.array([line.interpolate(s).coords[0] for s in distances])
        centre_levels=sample_surface(surface,xy,affine)
        slope=np.degrees(np.arctan(np.abs(np.gradient(centre_levels,distances))))
        for index,distance in enumerate(distances):
            point=Point(xy[index]);section=transect(line,distance,half_width)
            width,status=polygon_width(section,point,polygon)
            samples=np.array([section.interpolate(d+half_width).coords[0] for d in offsets])
            levels=sample_surface(surface,samples,affine)
            sid=f'{segment:02d}-{index:04d}'
            record=dict(station=sid,segment=segment,downstream_m=float(distance),
                        mapped_width_m=width,width_status=status,
                        original_centre_elevation_m=float(centre_levels[index]) if np.isfinite(centre_levels[index]) else None,
                        along_route_slope_deg=float(slope[index]) if np.isfinite(slope[index]) else None,
                        exceeds_slope_threshold=bool(slope[index]>maximum_slope) if np.isfinite(slope[index]) else None)
            stations.append(record)
            sections.append(dict(type='Feature',properties=record,geometry=mapping(transform(unproject,section))))
            cross_rows.extend((sid,float(offset),float(z) if np.isfinite(z) else None) for offset,z in zip(offsets,levels))
    output.mkdir(parents=True,exist_ok=True)
    with (output/'stations.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(stations[0]));writer.writeheader();writer.writerows(stations)
    with (output/'terrain-cross-sections.csv').open('w',newline='') as handle:
        writer=csv.writer(handle);writer.writerow(['station','offset_m','original_elevation_m']);writer.writerows(cross_rows)
    (output/'transects.geojson').write_text(json.dumps(dict(type='FeatureCollection',features=sections,attribution=snapshot['attribution']),indent=2)+'\n')
    valid=[r['mapped_width_m'] for r in stations if r['mapped_width_m'] is not None]
    report=dict(location=location,mode='diagnostic only; no terrain edits',river='Lower Hollyford',
                route_segments=len(lines),stations=len(stations),width_status_counts=dict(Counter(r['width_status'] for r in stations)),
                mapped_width_min_m=min(valid) if valid else None,mapped_width_median_m=float(np.median(valid)) if valid else None,
                mapped_width_max_m=max(valid) if valid else None,maximum_slope_deg=maximum_slope,
                settings=dict(sample_spacing_m=spacing,transect_half_width_m=half_width,tangent_window_m=25,cross_section_sample_m=5),
                inputs=dict(baseline_manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
                            lake_boundaries_sha256=hashlib.sha256(lake_path.read_bytes()).hexdigest(),
                            surface=str(surface_path),surface_sha256=hashlib.sha256(surface_path.read_bytes()).hexdigest(),
                            polygons_sha256=hashlib.sha256(SNAPSHOT.read_bytes()).hexdigest(),routes_sha256=hashlib.sha256(routes_path.read_bytes()).hexdigest(),
                            polygon_source=SERVICE),
                limitations=['Widths are mapped polygon cross-sections, not surveyed widths at a specified water level.',
                             'Outside-polygon, truncated and ambiguous sections remain unknown; no narrow-width fallback.',
                             'Disconnected map-edge segments are kept separate; islands are not bridged.',
                             'Cross-section elevations are original land DEM samples, not measured underwater beds.',
                             'Slope flags are diagnostics, not implemented excavation exclusions or source basins.',
                             'No tapering, bank detection, excavation decisions, water simulation or map export in this lab yet.'])
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    plot_review(output,lines,polygon,stations,cross_rows,sections,project)
    print(json.dumps(report,indent=2))
    return report


def plot_review(output,lines,polygon,stations,cross_rows,sections,project):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(13,7),layout='constrained')
    for boundary in line_parts(polygon.boundary):
        x,y=boundary.xy;axes[0].plot(x,y,color='#4894aa',linewidth=1)
    for line in lines:
        x,y=line.xy;axes[0].plot(x,y,color='#263b48',linewidth=1)
    for f in sections[::4]:
        x,y=transform(project,shape(f['geometry'])).xy;axes[0].plot(x,y,color='#df8231',alpha=.35,linewidth=.6)
    axes[0].set(aspect='equal',title='Lower Hollyford: mapped banks and transects',xlabel='Local easting (m)',ylabel='Local northing (m)')
    for segment in sorted({r['segment'] for r in stations}):
        rows=[r for r in stations if r['segment']==segment]
        axes[1].plot([r['downstream_m'] for r in rows],[r['mapped_width_m'] if r['mapped_width_m'] is not None else np.nan for r in rows],label=f'Segment {segment}')
    axes[1].set(title='Polygon-derived channel spans; gaps = unknown',xlabel='Distance along clipped segment (m)',ylabel='Mapped width (m)');axes[1].legend();axes[1].grid(alpha=.2)
    fig.savefig(output/'width-review.png',dpi=150);plt.close(fig)
    selected=[stations[int(i)]['station'] for i in np.linspace(0,len(stations)-1,min(4,len(stations)),dtype=int)]
    fig,axes=plt.subplots(len(selected),1,figsize=(10,3*len(selected)),layout='constrained',squeeze=False)
    for ax,sid in zip(axes.flat,selected):
        rows=[r for r in cross_rows if r[0]==sid]
        ax.plot([r[1] for r in rows],[r[2] for r in rows]);ax.axvline(0,color='gray',linestyle='--')
        ax.set(title=f'Original terrain — station {sid}',xlabel='Across route (m)',ylabel='Elevation (m)');ax.grid(alpha=.2)
    fig.savefig(output/'terrain-cross-sections.png',dpi=150);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['inspect','fetch-polygons'],nargs='?',default='inspect')
    parser.add_argument('--baseline-root',type=Path,default=terrain.BASE,help='Read-only checkout containing the existing Jamestown baseline/cache')
    parser.add_argument('--output',type=Path,default=terrain.BASE/'output/jamestown-river-lab')
    parser.add_argument('--sample-spacing',type=float,default=50.)
    parser.add_argument('--transect-half-width',type=float,default=500.)
    parser.add_argument('--maximum-slope',type=float,default=45.)
    args=parser.parse_args()
    if not all(np.isfinite(x) and x>0 for x in [args.sample_spacing,args.transect_half_width,args.maximum_slope]) or args.maximum_slope>=90:
        parser.error('Spacing/half-width must be finite and positive; slope must be between 0 and 90 degrees')
    if args.command=='fetch-polygons':fetch_polygons()
    else:inspect(args.baseline_root.resolve(),args.output.resolve(),args.sample_spacing,args.transect_half_width,args.maximum_slope)


if __name__=='__main__':main()
