"""Validate and plot the exported preset polygon trial, then refresh its ZIP."""
import argparse
import json,hashlib,shutil
from pathlib import Path
import numpy as np
from PIL import Image
from pyproj import Transformer
from shapely.geometry import shape
from shapely.ops import transform,unary_union
from nz_terrain import terrain,hydrology
import matplotlib.pyplot as plt


def review(baseline_root,slug='jamestown'):
    root=terrain.BASE;label=f'{slug}-polygon-rivers';out=root/'output'/label;base=baseline_root/'output'/slug
    location=json.loads((root/'locations.json').read_text())[slug]
    project=Transformer.from_crs(4326,terrain.projection(location['longitude'],location['latitude']),always_xy=True).transform
    features=json.loads((out/f'metadata/{slug}-polygons.geojson').read_text())['features']
    polys=[transform(project,shape(f['geometry'])) for f in features]
    opening_path=out/'metadata'/('jamestown-hollyford-openings.geojson' if slug=='jamestown' else f'{slug}-openings.geojson')
    opening_features=json.loads(opening_path.read_text())['features'] if opening_path.exists() else []
    opens=[transform(project,shape(f['geometry'])) for f in opening_features]
    allowed=unary_union(polys+[o.buffer(100) for o in opens])
    lakes=json.loads((out/'metadata/lake-polygons.geojson').read_text())['features']
    checks={};data={};changed={}
    for kind in terrain.WIDTHS:
        new=np.array(Image.open(out/f'heightmaps/{label}-{kind}.png'))
        old=np.array(Image.open(terrain.output_file(base,f'{slug}-{kind}.png')))
        assert new.dtype==np.uint16 and new.shape==(4096,4096)
        mask=np.zeros(new.shape,bool);r,c,m=hydrology.lake_window(allowed,terrain.grid(location,terrain.WIDTHS[kind]),4096);mask[r:r+m.shape[0],c:c+m.shape[1]]=m
        changed[kind]=new!=old
        if kind=='world':mask[1536:2560,1536:2560]|=changed['playable'].reshape(1024,4,1024,4).any(axis=(1,3))
        assert not (changed[kind]&~mask).any()
        if kind=='playable':assert np.all(new<=old)
        checks[kind]=dict(changed_pixels=int(changed[kind].sum()),outside_water_and_nested_sync_changed_pixels=int((changed[kind]&~mask).sum()),max_additional_lowering_m=float((old.astype(float)-new).max()*terrain.STEP))
        data[kind]=new
        terrain.preview(new,location,kind,out/f'previews/{kind}-preview.jpg',lake_mask=hydrology.feature_mask(terrain,location,kind,lakes))
    affected=changed['playable'].reshape(1024,4,1024,4).any(axis=(1,3))
    means=data['playable'].astype(float).reshape(1024,4,1024,4).mean(axis=(1,3))
    err=float(np.abs(data['world'][1536:2560,1536:2560][affected]-means[affected]).max()*terrain.STEP) if affected.any() else 0.
    assert err<=terrain.STEP/2+1e-10
    checks['nested_grid_max_error_m']=float(err)
    checks['lakes']={}
    for f in lakes:
        if any(n in (f['properties'].get('name') or '').lower() for n in ['mckerrow','waiuna']):
            geom=transform(project,shape(f['geometry']));r,c,m=hydrology.lake_window(geom,terrain.grid(location,14336),4096)
            checks['lakes'][f['properties']['name']]={'playable_pixels':int(m.sum()),'lowest_bed_elevation_m':float(data['playable'][r:r+m.shape[0],c:c+m.shape[1]][m].min()*terrain.STEP+terrain.HEIGHT_MIN)}
    report_path=out/'metadata/river-report.json';report=json.loads(report_path.read_text());report['validation']=checks
    report['changed_pixels']={k:v['changed_pixels'] for k,v in checks.items() if k in terrain.WIDTHS}
    report['baseline_heightmap_sha256']={k:hashlib.sha256(terrain.output_file(base,f'{slug}-{k}.png').read_bytes()).hexdigest() for k in terrain.WIDTHS}
    settings=json.loads((root/'water-configs/river-preferred.json').read_text());report['settings']=settings
    report_path.write_text(json.dumps(report,indent=2)+'\n');shutil.copy(root/'water-configs/river-preferred.json',out/'metadata/river-preferred.json')
    if slug=='jamestown':
        # Show actual exported elevations at both openings, on identical colour scales.
        fig,axes=plt.subplots(2,2,figsize=(10,10),layout='constrained')
        aff=terrain.grid(location,14336)
        old=np.array(Image.open(terrain.output_file(base,'jamestown-playable.png')))
        for row,opening in enumerate(opens):
            pt=opening.interpolate(.5,normalized=True);c,r=(~aff)*(pt.x,pt.y);r=int(r);c=int(c);n=90
            r0=max(0,r-n);r1=min(4096,r+n);c0=max(0,c-n);c1=min(4096,c+n)
            for col,(a,label) in enumerate([(old,'Existing lake/sea baseline'),(data['playable'],'Polygon river trial')]):
                z=a[r0:r1,c0:c1].astype(float)*terrain.STEP+terrain.HEIGHT_MIN
                im=axes[row,col].imshow(z,cmap='terrain',vmin=-12,vmax=12,extent=(0,(c1-c0)*3.5,(r1-r0)*3.5,0))
                axes[row,col].set(title=f'{["Sea mouth","Lake McKerrow connection"][row]}\n{label}',xlabel='Map east (m)',ylabel='Map south (m)')
        fig.colorbar(im,ax=axes,label='Exported bed / land elevation (m)',shrink=.6)
        fig.suptitle('Hollyford connections in the exported playable heightmap')
        fig.savefig(out/'previews/hollyford-connections.png',dpi=160);plt.close(fig)
    terrain.package_outputs(out,label)
    print(json.dumps(checks,indent=2))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-root",type=Path,required=True)
    parser.add_argument('--location',default='jamestown')
    args=parser.parse_args()
    review(args.baseline_root.resolve(),args.location)


if __name__=="__main__":main()
