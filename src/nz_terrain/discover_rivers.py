"""Snapshot public LINZ river centrelines for the configured playable footprint.

Run before nz-terrain --with-rivers. Existing OSM supplements are retained;
wide rivers missing from LINZ's line layer still need explicit supplements.
"""
import argparse
import json
from urllib.parse import urlencode
from shapely.geometry import shape, Polygon
from . import terrain

SERVICE = 'https://services.arcgis.com/xdsHIIxuCWByZiCB/ArcGIS/rest/services/LINZ_NZ_River_Name_Lines_Pilot/FeatureServer/0'


def discover(slug):
    location = json.loads((terrain.BASE/'locations.json').read_text())[slug]
    boundary = Polygon(terrain.outline(location, terrain.WIDTHS['playable']))
    bbox = boundary.bounds
    features = []
    seen = set()
    offset = 0
    while True:
        query = dict(where='1=1', geometry=','.join(map(str,bbox)), geometryType='esriGeometryEnvelope',
                     inSR=4326, outSR=4326, spatialRel='esriSpatialRelIntersects', outFields='*',
                     f='geojson', resultOffset=offset, resultRecordCount=1000, orderByFields='OBJECTID')
        page = terrain.fetch(SERVICE+'/query?'+urlencode(query))
        if 'error' in page: raise RuntimeError(page['error'])
        rows = page.get('features', [])
        for f in rows:
            key = f['properties']['OBJECTID']
            if key in seen: raise RuntimeError('Repeated page from river service')
            seen.add(key)
            if f['geometry'] and shape(f['geometry']).intersects(boundary):
                f['properties']['source'] = 'LINZ River Name Lines Pilot, CC BY 4.0'
                features.append(f)
        offset += len(rows)
        more=page.get('exceededTransferLimit') or page.get('properties',{}).get('exceededTransferLimit')
        if not rows and more:raise RuntimeError('River service pagination stalled')
        if not rows or (len(rows)<1000 and not more): break
    path = terrain.BASE/'river-configs'/f'{slug}-network.geojson'
    supplements = []
    if path.exists():
        supplements = [f for f in json.loads(path.read_text())['features']
                       if 'OpenStreetMap' in f['properties'].get('source','') and shape(f['geometry']).intersects(boundary)]
    names = {f['properties'].get('name') for f in supplements if f['properties'].get('name')}
    if 'Tūtoko River' in names:names.add('Tutoko River')
    features = [f for f in features if f['properties'].get('name') not in names] + supplements
    data = dict(type='FeatureCollection', source=SERVICE, query_bbox=list(bbox), location=location,
                coverage_note='LINZ line features omit some wide river channels and do not cover all offshore islands. Widths are inferred, not measured.', features=features)
    path.write_text(json.dumps(data,ensure_ascii=False)+'\n')
    configpath = terrain.BASE/'river-configs'/f'{slug}.json'
    if not configpath.exists() and features:
        config = json.loads((terrain.BASE/'river-configs/jamestown.json').read_text())
        config.update(location=slug,routes=path.name,named_overrides={})
        configpath.write_text(json.dumps(config,indent=2)+'\n')
    print(slug, len(features), 'features;', sorted({f['properties'].get('name') or '(unnamed)' for f in features}),flush=True)
    return data


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--location',action='append')
    parser.add_argument('--all',action='store_true')
    args=parser.parse_args()
    if args.all==bool(args.location):parser.error('Choose --all or --location')
    for slug in (json.loads((terrain.BASE/'locations.json').read_text()) if args.all else args.location):
        discover(slug)


if __name__ == '__main__':
    main()
