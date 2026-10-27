"""Inspect Hollyford end caps and demonstrate explicit bank exclusions.

Uses the existing river snapshot and OSM route, with a read-only lake baseline.
No general-purpose automatic opening classifier or production integration yet.
"""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer
from shapely.geometry import LineString, MultiLineString, shape, mapping
from shapely.ops import transform, unary_union
from water_carver import BedParameters, generate_bed
from plot_samples import ROOT, grid_for
from plot_two_pass import profile, smooth_creases


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-root',type=Path,required=True)
    args=parser.parse_args()
    out=ROOT/'output/jamestown-river-lab/openings';out.mkdir(parents=True,exist_ok=True)
    source=ROOT/'river-configs/jamestown-hollyford-polygons.geojson'
    project=Transformer.from_crs(4326,2193,always_xy=True).transform
    unproject=Transformer.from_crs(2193,4326,always_xy=True).transform
    source_poly=shape(json.loads(source.read_text())['features'][0]['geometry'])
    water=transform(project,source_poly)
    route_source=ROOT/'river-configs/jamestown-network.geojson'
    route=transform(project,shape(next(f for f in json.loads(route_source.read_text())['features'] if 'hollyford' in f['properties'].get('name','').lower())['geometry']))
    lakepath=args.baseline_root/'output/jamestown/metadata/lake-polygons.geojson'
    lakefeatures=json.loads(lakepath.read_text())['features']
    lake=unary_union([transform(project,shape(f['geometry'])) for f in lakefeatures if 'mckerrow' in (f['properties'].get('name') or '').lower()])
    if lake.is_empty:raise ValueError('Missing Lake McKerrow baseline')
    coordinates=list(water.exterior.coords);original_coords=list(source_poly.exterior.coords)
    segments=[LineString([a,b]) for a,b in zip(coordinates[:-1],coordinates[1:])]
    candidates=[]
    for i,edge in enumerate(segments):
        intersection=edge.intersection(route)
        if intersection.is_empty or intersection.geom_type!='Point':continue
        station=route.project(intersection)
        before=route.interpolate(max(0,station-30));after=route.interpolate(min(route.length,station+30))
        a=np.array(after.coords[0])-np.array(before.coords[0]);b=np.diff(np.array(edge.coords),axis=0)[0]
        angle=np.degrees(np.arccos(np.clip(abs(a@b)/(np.linalg.norm(a)*np.linalg.norm(b)),0,1)))
        lake_shared=edge.intersection(lake.boundary).length
        label='Lake McKerrow connection' if lake_shared>edge.length*.95 else 'Sea-mouth candidate'
        props=dict(name=label,exterior_segment_index=i,length_m=edge.length,route_crossing_angle_deg=float(angle),lake_boundary_shared_m=lake_shared,route_distance_to_downstream_end_m=route.length-station,status='diagnostic candidate; not used by production generation',evidence='Exact shared boundary with Lake McKerrow' if lake_shared>edge.length*.95 else 'Exterior crossing near downstream OSM route endpoint; coastal location. Not verified against a separate sea polygon.')
        candidates.append((props,edge))
    if len(candidates)!=2:raise ValueError('Expected two exterior crossings; review updated inputs rather than guessing')
    excluded={p['exterior_segment_index'] for p,_ in candidates}
    # Original remaining segments plus every island shore; never buffer/remove banks around a crossing.
    banks=MultiLineString([list(g.coords) for i,g in enumerate(segments) if i not in excluded]+[list(r.coords) for r in water.interiors])
    np.testing.assert_allclose(water.boundary.length-banks.length,sum(g.length for _,g in candidates),atol=1e-6)
    for hole in water.interiors:
        np.testing.assert_allclose(banks.intersection(LineString(hole.coords)).length, hole.length, atol=1e-6)
    features=[dict(type='Feature',properties=p,geometry=mapping(LineString([original_coords[p['exterior_segment_index']],original_coords[p['exterior_segment_index']+1]]))) for p,_ in candidates]
    snap=dict(type='FeatureCollection',features=features,source=str(source.relative_to(ROOT)),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),route_sha256=hashlib.sha256(route_source.read_bytes()).hexdigest(),lake_baseline_sha256=hashlib.sha256(lakepath.read_bytes()).hexdigest(),attribution='LINZ polygon CC BY 4.0; crossing evidence OpenStreetMap contributors ODbL 1.0',scope='Two Hollyford diagnostic opening candidates, not an automatic rule for other rivers.')
    (out/'opening-candidates.geojson').write_text(json.dumps(snap,indent=2)+'\n')
    (out/'retained-banks.geojson').write_text(json.dumps(dict(type='FeatureCollection',features=[dict(type='Feature',properties={'role':'bank','status':'diagnostic'},geometry=mapping(transform(unproject,banks)))]))+'\n')
    settings=json.loads((ROOT/'water-configs/river-preferred.json').read_text())
    fig,axes=plt.subplots(2,3,figsize=(15,10),layout='constrained')
    for row,(props,edge) in enumerate(candidates):
        pt=edge.interpolate(.5,normalized=True);half=260
        grid=grid_for((pt.x-half,pt.y-half,pt.x+half,pt.y+half),1)
        ax=axes[row,0]
        for ring in [water.exterior,*water.interiors]:
            x,y=ring.xy;ax.plot(np.asarray(x)-pt.x,np.asarray(y)-pt.y,color='#6b8590',lw=1)
        if props['lake_boundary_shared_m']>0:
            parts=list(lake.geoms) if hasattr(lake,'geoms') else [lake]
            for l in parts:
                x,y=l.exterior.xy;ax.plot(np.asarray(x)-pt.x,np.asarray(y)-pt.y,color='#409fc9',lw=1,label='Lake boundary')
        x,y=route.xy;ax.plot(np.asarray(x)-pt.x,np.asarray(y)-pt.y,color='#799644',lw=1,label='Existing OSM route')
        x,y=edge.xy;ax.plot(np.asarray(x)-pt.x,np.asarray(y)-pt.y,color='#cf553b',lw=3,label='Excluded cap')
        ax.set(title=f"{props['name']}\n{edge.length:.1f} m boundary segment",xlim=(-half,half),ylim=(-half,half),aspect='equal',xlabel='East of opening midpoint (m)',ylabel='North of opening midpoint (m)');ax.legend(fontsize=8,loc='lower left')
        for col,(bank_geom,title) in enumerate([(None,'All boundary edges treated as banks'),(banks,'Opening excluded from bank distance')],1):
            bed=generate_bed(water,grid,water_level_m=0,parameters=BedParameters(10000,1,'linear'),banks=bank_geom)
            z=profile(bed.depth_m,settings['shore_transition_m'],cap=settings['max_depth_m'],shore_gradient=settings['shore_gradient'],deep_gradient=settings['interior_gradient'])
            s=settings['smoothing'];z,_=smooth_creases(z,bed.depth_m,1,s['scale_m'],s['shore_lock_m'],s['shore_blend_m'])
            im=axes[row,col].imshow(z,extent=(-half,half,-half,half),origin='upper',cmap='YlGnBu',vmin=0,vmax=10,interpolation='nearest')
            x,y=edge.xy;axes[row,col].plot(np.asarray(x)-pt.x,np.asarray(y)-pt.y,color='#cf553b',lw=1,ls='--')
            axes[row,col].set(title=title,xlabel='East of opening midpoint (m)',ylabel='North of opening midpoint (m)')
    fig.colorbar(im,ax=list(axes[:,1:].flat),shrink=.55,label='Synthetic depth (m)')
    fig.suptitle('Hollyford end caps / remove artificial banks, retain the water footprint\nPreferred river profile: 10 m cap · 2 m transition · 1 m smoothing · 100% nominal bank slope',fontsize=17)
    fig.text(.5,.005,'Diagnostic only: flat assumed water surface at 0 m; adjacent lake/sea beds are not joined or modified. LINZ CC BY 4.0; OSM ODbL.',ha='center',fontsize=9)
    fig.set_constrained_layout_pads(h_pad=.4,w_pad=.1)
    fig.savefig(out/'opening-comparison.png',dpi=170);plt.close(fig)
    report=dict(exterior_opening_candidates=[p for p,_ in candidates],island_boundary_route_crossings=sum(len(g.geoms) if hasattr(g,'geoms') else int(not g.is_empty) for g in [route.intersection(LineString(r.coords)) for r in water.interiors]),remaining_bank_length_m=banks.length,settings=settings)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
