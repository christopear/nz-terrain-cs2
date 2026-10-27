"""Save LINZ river-area polygons intersecting a preset's playable envelope."""
import argparse
from datetime import datetime, timezone
import json
from urllib.parse import urlencode
from urllib.request import urlopen
from .lab import SERVICE
from .. import terrain


def fetch(slug):
    location=json.loads((terrain.BASE/'locations.json').read_text())[slug]
    params=dict(f='geojson',where='1=1',geometry=','.join(map(str,terrain.bounds(terrain.outline(location,14336)))),
                geometryType='esriGeometryEnvelope',spatialRel='esriSpatialRelIntersects',inSR=4326,outSR=4326,
                outFields='*',resultRecordCount=1000,orderByFields='OBJECTID')
    features=[];seen=set();urls=[]
    while True:
        params['resultOffset']=len(features)
        url=SERVICE+'/query?'+urlencode(params);urls.append(url)
        with urlopen(url,timeout=60) as response:page=json.load(response)
        if 'features' not in page:raise ValueError(f'Invalid polygon response: {page}')
        batch=page['features']
        for feature in batch:
            key=feature['properties']['OBJECTID']
            if key in seen:raise ValueError('Repeated polygon page')
            seen.add(key);features.append(feature)
        more=page.get('exceededTransferLimit') or page.get('properties',{}).get('exceededTransferLimit')
        if not batch and more:raise ValueError('Pagination stalled')
        if not more and len(batch)<1000:break
    result=dict(type='FeatureCollection',features=features,location=location,source=SERVICE,query_urls=urls,
                retrieved_utc=datetime.now(timezone.utc).isoformat(),coverage='Playable-map geographic envelope; full uncut source polygons retained',
                attribution='LINZ River Name Polygons Pilot, CC BY 4.0. Mapped water areas, not surveyed bathymetry.')
    destination=terrain.BASE/'river-configs'/f'{slug}-polygons.geojson'
    temporary=destination.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(destination)
    print(slug,len(features),[(f['properties'].get('OBJECTID'),f['properties'].get('name')) for f in features],flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--location',action='append',required=True)
    for slug in parser.parse_args().location:fetch(slug)


if __name__=='__main__':main()
