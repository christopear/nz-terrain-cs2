"""Experimental, cut-only river beds. No water simulation or game water sources."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
os.environ.setdefault('MPLCONFIGDIR', '/tmp/jamestown-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from PIL import Image
from pyproj import Transformer
from scipy.ndimage import gaussian_filter1d, median_filter
from shapely.geometry import shape, LineString, Point, Polygon, mapping
from shapely.ops import transform, substring
from . import terrain
from .river.core import solve_profile, sample, carve as carve_grid




def load_routes(config, location):
    crs = terrain.projection(location['longitude'],location['latitude'])
    to_local = Transformer.from_crs(4326,crs,always_xy=True).transform
    footprint = Polygon([terrain.grid(location,14336)*p for p in [(0,0),(4096,0),(4096,4096),(0,4096)]])
    features = json.loads((terrain.BASE/'river-configs'/config['routes']).read_text())['features']
    lines = {}
    for f in features:
        line = transform(to_local,shape(f['geometry'])).intersection(footprint)
        if line.geom_type != 'LineString':
            raise ValueError(f"Route {f['properties']['id']} needs explicit selection after clipping: {line.geom_type}")
        trim=config['rivers'][f['properties']['id']].get('trim_upstream_m',0)
        if trim:line=substring(line,trim,line.length)
        lines[f['properties']['id']] = line
    main = lines['takaka']
    junctions = {}
    for name,line in lines.items():
        if name=='takaka':continue
        end = Point(line.coords[-1]); s = main.project(end); p=main.interpolate(s)
        if p.distance(end)>5:raise ValueError(f'{name} disconnected from main river by {p.distance(end):.1f} m')
        lines[name]=LineString(list(line.coords)[:-1]+[p.coords[0]])
        junctions[name]=s
    return lines,junctions




def analyse(config,location,z):
    lines,junctions=load_routes(config,location)
    profiles={}
    affine=terrain.grid(location,14336)
    for name,line in lines.items():
        s=np.linspace(0,line.length,int(np.ceil(line.length/config['sample_spacing_m']))+1)
        # Include actual vertices and main-river junctions to avoid cutting corners.
        s=np.unique(np.r_[s,[line.project(Point(p)) for p in line.coords],list(junctions.values()) if name=='takaka' else []])
        s=s[np.r_[True,np.diff(s)>1e-5]]
        xy=np.array([line.interpolate(d).coords[0] for d in s])
        centre=sample(z,xy,affine)
        # Sample a small corridor so mapped-centreline error doesn't turn an
        # individual gravel bar into a spurious dam. Never search across valleys.
        offsets=np.array([(x,y) for x in [-7,0,7] for y in [-7,0,7]])
        low=np.percentile(np.array([sample(z,xy+o,affine) for o in offsets]),25,axis=0)
        # Resample uniformly for metre-based filtering (route vertex spacing varies).
        uniform=np.arange(0,line.length+5,5.)
        ref=np.interp(s,uniform,gaussian_filter1d(median_filter(np.interp(uniform,s,low),size=5),config['surface_smoothing_m']/5))
        ref=np.maximum(ref,0)
        profiles[name]={'s':s,'xy':xy,'centre':centre,'reference':ref}
    for name in ['takaka']+[n for n in profiles if n!='takaka']:
        p=profiles[name]
        outlet=config['outlet_bed_m'] if name=='takaka' else float(np.interp(junctions[name],profiles['takaka']['s'],profiles['takaka']['bed']))
        p['bed']=solve_profile(p['reference'],p['s'],config['rivers'][name]['depth_m'],config['max_cut_m'],config['preferred_bed_slope'],outlet)
    return profiles,lines,junctions


def carve(z,profiles,config,location):
    """Compatibility adapter for the existing map-generation pipeline."""
    return carve_grid(z, profiles, config, terrain.grid(location, 14336))


def plots(out,z,carved,profiles,location):
    affine=terrain.grid(location,14336)
    fig,ax=plt.subplots(figsize=(11,10),layout='constrained')
    small=z[::4,::4];dy,dx=np.gradient(small,14)
    shade=np.clip(.65+(-dx+dy)*.2,.15,1)
    ax.imshow(shade,cmap='gray',vmin=0,vmax=1,extent=(0,14336,14336,0))
    water=small<0
    if (out/'lake-polygons.geojson').exists():
        from . import hydrology
        water|=hydrology.feature_mask(terrain,location,'playable',json.loads((out/'lake-polygons.geojson').read_text())['features'])[::4,::4]
    rgba=np.zeros((*water.shape,4));rgba[:,:,:3]=[.15,.55,.8];rgba[:,:,3]=water*.35
    ax.imshow(rgba,extent=(0,14336,14336,0))
    colors={};label_colors={}
    seen=set()
    for name,p in profiles.items():
        label=p.get('label',name);legend=label if label not in seen and len(seen)<12 else '_nolegend_';seen.add(label)
        colour=label_colors.setdefault(label,plt.get_cmap('tab10')(len(label_colors)%10))
        c,r=(~affine)*(p['xy'][:,0],p['xy'][:,1]);line,=ax.plot(c*3.5,r*3.5,lw=1.2,label=legend,color=colour);colors[name]=line.get_color()
        marker=min(5,len(c)-1)
        ax.scatter(c[marker]*3.5,r[marker]*3.5,s=15 if len(profiles)>10 else 50,color=colors[name],marker='v')
    ax.set(title=location['name']+' river candidate — routes and upstream source markers',xlabel='Metres from left map edge',ylabel='Metres from top map edge');ax.legend(fontsize=8);fig.savefig(out/'river-overview.png',dpi=150);plt.close(fig)
    selected=sorted(profiles.items(),key=lambda item:item[1]['s'][-1],reverse=True)[:6]
    fig,axes=plt.subplots(len(selected),1,figsize=(12,3*len(selected)),layout='constrained',squeeze=False)
    for ax,(name,p) in zip(axes.flat,selected):
        after=sample(carved,p['xy'],affine)
        ax.plot(p['s']/1000,p['centre'],color='.65',lw=.8,label='Baseline at mapped centreline')
        ax.plot(p['s']/1000,p['reference'],color='steelblue',lw=1,label='Smoothed surface proxy (not water simulation)')
        ax.plot(p['s']/1000,p['bed'],color='darkorange',lw=1,label='Designed bed')
        ax.plot(p['s']/1000,after,color='black',lw=.6,label='Exported terrain along route')
        ax.set(title=p.get('label',name)+' — '+name,xlabel='Downstream distance (km)',ylabel='Elevation (m)');ax.grid(alpha=.2)
    axes.flat[0].legend(fontsize=8);fig.savefig(out/'river-profiles.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for ax,(name,p) in zip(axes.flat,selected[:4]):
        idx=len(p['s'])//2;tangent=p['xy'][idx+1]-p['xy'][idx-1];normal=np.array([-tangent[1],tangent[0]])/np.linalg.norm(tangent)
        offset=np.linspace(-50,50,201);xy=p['xy'][idx]+offset[:,None]*normal
        ax.plot(offset,sample(z,xy,affine),label='Baseline');ax.plot(offset,sample(carved,xy,affine),label='Candidate')
        ax.set(title=p.get('label',name)+' — mid-reach cross-section',xlabel='Across river (m)',ylabel='Elevation (m)');ax.grid(alpha=.2)
    axes.flat[0].legend();fig.savefig(out/'river-cross-sections.png',dpi=150);plt.close(fig)


def build(slug,options=None):
    config=json.loads((terrain.BASE/'river-configs'/f'{slug}.json').read_text())
    if options and config.get('mode')!='network':
        raise ValueError('River width/depth CLI defaults apply to network mode; edit the per-river config for Takaka')
    config.update(options or {})
    location=json.loads((terrain.BASE/'locations.json').read_text())[slug]
    baseline=terrain.BASE/'output'/slug
    manifest=json.loads((terrain.output_file(baseline,'manifest.json')).read_text())
    if manifest['location']!=location:raise ValueError('Rebuild baseline for current location first')
    original=np.array(Image.open(terrain.output_file(baseline,f'{slug}-playable.png')))
    if 'source_surface_cache' not in manifest:
        raise ValueError('Rebuild baseline first: original source surface is required to avoid stacked excavation')
    with rasterio.open(terrain.BASE/manifest['source_surface_cache']/'surface-playable.tif') as src:
        source=terrain.encode(src.read(1))
    z=source.astype('float32')*terrain.STEP+terrain.HEIGHT_MIN
    if config.get('mode')=='network':
        from . import river_network
        config['rivers']={}
        profiles,lines,junctions=river_network.analyse(config,location,z,baseline)
    else:profiles,lines,junctions=analyse(config,location,z)
    carved,mask=carve(z,profiles,config,location)
    river_target=terrain.encode(carved)
    encoded=original.copy()
    encoded[mask]=np.minimum(original[mask],river_target[mask])
    assert np.array_equal(encoded[~mask],original[~mask])
    actual=encoded.astype('float32')*terrain.STEP+terrain.HEIGHT_MIN
    river_cut=z-(river_target.astype('float32')*terrain.STEP+terrain.HEIGHT_MIN)
    cut=(original.astype('float32')-encoded)*terrain.STEP
    allowed_cut=max(d.get('max_cut_m',config['max_cut_m']) for d in config['rivers'].values())
    if river_cut.max()>allowed_cut+terrain.STEP:
        raise ValueError(f'River target excavation {river_cut.max():.2f} m exceeds cut limit; inspect route alignment')
    world=np.array(Image.open(terrain.output_file(baseline,f'{slug}-world.png')))
    changed=encoded!=original
    affected=changed.reshape(1024,4,1024,4).any(axis=(1,3))
    means=encoded.astype('float64').reshape(1024,4,1024,4).mean(axis=(1,3))
    centre=world[1536:2560,1536:2560];centre[affected]=np.rint(means[affected]).astype('uint16')
    out=terrain.BASE/'output'/f'.{slug}-rivers-building'
    if out.exists():shutil.rmtree(out)
    out.mkdir()
    for filename in ['water-report.json','lake-polygons.geojson','water-profiles.png']:
        source_file=terrain.output_file(baseline,filename)
        if source_file.exists():shutil.copy(source_file,out/filename)
    stats={};features=[]
    to_geo=Transformer.from_crs(terrain.projection(location['longitude'],location['latitude']),4326,always_xy=True)
    for name,p in profiles.items():
        slopes=-np.diff(p['bed'])/np.diff(p['s'])
        xy=p['xy'][min(5,len(p['xy'])-1)];lon,lat=to_geo.transform(*xy)
        measured=sample(actual,p['xy'],terrain.grid(location,14336))
        dense_s=np.linspace(0,lines[name].length,int(np.ceil(lines[name].length/2))+1)
        dense_xy=np.array([lines[name].interpolate(s).coords[0] for s in dense_s])
        dense_bed=sample(actual,dense_xy,terrain.grid(location,14336))
        stats[name]={'length_m':float(p['s'][-1]),'upstream_bed_m':float(p['bed'][0]),'downstream_bed_m':float(p['bed'][-1]),'preferred_slope_shortfall_length_m':float(np.diff(p['s'])[slopes<config['preferred_bed_slope']-1e-6].sum()),'maximum_exported_uphill_step_m':float(max(0,np.diff(measured).max())),'suggested_source':{'longitude':lon,'latitude':lat,'surface_proxy_m':float(p['reference'][min(5,len(p['xy'])-1)]),'note':'Manual source placement candidate, not a calibrated inflow. Border reaches need interior sources; connected tributaries and lake outlets do not each require an extra source.'}}
        stats[name]['maximum_uphill_step_at_2m_sampling_m']=float(max(0,np.diff(dense_bed).max()))
        stats[name]['maximum_above_designed_bed_at_2m_sampling_m']=float(max(0,np.max(dense_bed-np.interp(dense_s,p['s'],p['bed']))))
        stats[name]['cut_limit_pixels_in_map']=int(np.count_nonzero(cut>=config['max_cut_m']-terrain.STEP)) if name=='takaka' else None
        if stats[name]['maximum_uphill_step_at_2m_sampling_m']>.2 and config.get('mode')!='network':
            raise ValueError(f'{name}: exported channel has an uphill obstruction; inspect profile before publishing')
        features.append({'type':'Feature','properties':{'river':name,**config['rivers'][name]},'geometry':mapping(transform(to_geo.transform,lines[name]))})
        np.savetxt(out/f'{name}-profile.csv',np.c_[p['s'],p['centre'],p['reference'],p['bed'],measured],delimiter=',',header='downstream_m,baseline_m,surface_proxy_m,designed_bed_m,exported_bed_m',comments='')
    for name,a in [('playable',encoded),('world',world)]:
        path=out/f'{slug}-rivers-{name}.png';Image.fromarray(a).save(path)
        check=np.array(Image.open(path));assert check.dtype==np.uint16 and check.shape==(4096,4096) and np.array_equal(a,check)
        terrain.preview(a,location,name,out/f'{name}-preview.jpg')
    plots(out,original.astype('float32')*terrain.STEP+terrain.HEIGHT_MIN,actual,profiles,location)
    report={'location':location,'config':config,'baseline_sha256':hashlib.sha256((terrain.output_file(baseline,f'{slug}-playable.png')).read_bytes()).hexdigest(),'changed_pixels':int(changed.sum()),'max_cut_m':float(cut.max()),'excavation_m3':float(cut.sum()*3.5**2),'outside_corridor_changed_pixels':0,'rivers':stats,'junctions_main_distance_m':junctions,'limitations':['Experimental synthetic beds, not surveyed bathymetry.','No in-game or hydraulic simulation validation; terrain does not place water sources.','Mapped river centrelines may differ from current channels.','Preferred bed slope is soft, not a guarantee of water-surface slope.','Lake/lagoon beds inherited from the baseline; only configured/selected river reaches carved. See selection decisions and exclusions.','Coastline may change inside the explicit river corridor; preserved outside it.']}
    report['composition']='minimum of baseline water bed and independent river target, inside river corridor only'
    report['network_selection']=config.get('network_report')
    report['raster_checks']={'uphill_steps_over_0_2m':[name for name,stats in stats.items() if stats['maximum_uphill_step_at_2m_sampling_m']>.2], 'note':'Small raster/confluence depressions may produce uphill steps despite monotonic design profiles. Review flagged reaches in-game.'}
    (out/'river-report.json').write_text(json.dumps(report,indent=2))
    (out/'river-routes.geojson').write_text(json.dumps({'type':'FeatureCollection','features':features},indent=2))
    (out/'IMPORT.txt').write_text((baseline/'IMPORT.txt').read_text().replace(f'{slug}-playable.png',f'{slug}-rivers-playable.png').replace(f'{slug}-world.png',f'{slug}-rivers-world.png')+'\nRIVER EXPERIMENT\nRiver targets use the original land surface; the lower absolute target wins where rivers meet sea, lakes or lagoons.\nUse BOTH rivers PNGs together. The baseline remains in the sibling location folder.\nSee river-report.json, route overview, profiles and cross-sections.\nCreate upstream water sources manually, initially at the marked interior points.\nStart low and adjust inflow until water stays roughly 1.5–3 m above the carved bed\nwithout overtopping banks. No universal source-strength setting is assumed.\nAllow the simulation to settle; check each confluence, flat reach and estuary.\nIf shallow water stalls, tune supply/width before increasing slope.\nStop and compare baseline if sustained flooding or a channel discontinuity occurs.\nBed depths in config are excavation below a surface proxy, NOT guaranteed water depth.\nOnly configured/selected river reaches are included. Consult river-report.json for the route inventory and exclusions.\n\nRiver routes: © OpenStreetMap contributors, ODbL 1.0, via WaterwayMap.org.\nhttps://www.openstreetmap.org/copyright\nhttps://opendatacommons.org/licenses/odbl/1-0/\nhttps://waterwaymap.org/\nRiver route database and adaptations are available under ODbL 1.0.\n')
    shutil.copy(terrain.output_file(baseline,'manifest.json'),out/'baseline-manifest.json')
    shutil.copy(terrain.BASE/'river-configs/SOURCES.md',out/'river-sources.md')
    with (out/'IMPORT.txt').open('a') as attribution:
        attribution.write('\nAdditional river geometry: LINZ River Name Lines Pilot, CC BY 4.0.\nSee metadata/river-sources.md and the supplied route database for source details.\n')
    (out/'river-report.json').write_text(json.dumps(report,indent=2))
    terrain.package_outputs(out,f'{slug}-rivers')
    destination=terrain.BASE/'output'/f'{slug}-rivers'
    if destination.exists():shutil.rmtree(destination)
    out.replace(destination)
    print(json.dumps(report,indent=2))

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--location',required=True)
    parser.add_argument('--minimum-river-width',type=float,help='Estimated width cutoff in metres, network mode; trims headwaters')
    parser.add_argument('--default-river-depth',type=float,help='Base excavation depth in metres, network mode')
    args=parser.parse_args();options={}
    for value,key in [(args.minimum_river_width,'minimum_river_width_m'),(args.default_river_depth,'default_depth_m')]:
        if value is not None:
            if not np.isfinite(value) or value<=0:parser.error('Width and depth must be finite and positive')
            options[key]=value
    build(args.location,options)


if __name__ == '__main__':
    main()
