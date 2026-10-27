"""Keep only the user-reviewed Awarua and McKenzie polygon gap connections."""
import json
from nz_terrain import terrain

REVIEWED_PAIRS={(180,2238),(2238,2239),(2239,4562),(2240,4563),(1108,4563)}


def main():
    source=terrain.BASE/'output/jamestown-river-lab/network-trial'
    destination=terrain.BASE/'output/jamestown-river-lab/reviewed-connectors'
    report=json.loads((source/'connection-report.json').read_text())
    accepted=[];rejected=[]
    for item in report['generated']:
        if item['kind']=='polygon_connector' and tuple(sorted(item['polygon_ids'])) in REVIEWED_PAIRS:
            accepted.append(item)
        else:rejected.append(dict(item,status='excluded',reason='Excluded by user: retain only reviewed Awarua and McKenzie polygon connectors'))
    assert len(accepted)==7, 'Expected the two reviewed Awarua and five McKenzie channels'
    ids={item['id'] for item in accepted}
    destination.mkdir(parents=True,exist_ok=True)
    for name in ('synthetic-connectors.geojson','connector-traces.geojson'):
        data=json.loads((source/name).read_text());data['features']=[f for f in data['features'] if f['properties']['id'] in ids]
        assert len(data['features'])==7
        (destination/name).write_text(json.dumps(data)+'\n')
    report['generated']=accepted;report['skipped'].extend(rejected)
    report['scope']='Reviewed Awarua and McKenzie polygon-to-polygon connectors only; no inferred tributaries or new sea/lake links.'
    (destination/'connection-report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(destination, '—',len(accepted),'reviewed connectors')

if __name__=='__main__':main()
