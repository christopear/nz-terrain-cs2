"""Compare Waiuna shoreline transitions and distance smoothing; leaves map packages/settings unchanged."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from water_carver.river import profile
import numpy as np
from pyproj import Transformer
from scipy.ndimage import distance_transform_edt, gaussian_filter
from shapely.geometry import shape
from shapely.ops import transform
from nz_terrain import terrain, hydrology


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--short-transitions',action='store_true')
    args=parser.parse_args()
    root=terrain.BASE
    folder='short-shoreline-trials' if args.short_transitions else 'shoreline-trials'
    out=root/'output/jamestown-river-lab/waiuna'/folder;out.mkdir(parents=True,exist_ok=True)
    source=root/'output/jamestown-polygon-rivers/metadata/lake-polygons.geojson'
    feature=next(f for f in json.loads(source.read_text())['features'] if f['properties'].get('name')=='Waiuna Lagoon')
    location=json.loads((root/'locations.json').read_text())['jamestown']
    project=Transformer.from_crs(4326,terrain.projection(location['longitude'],location['latitude']),always_xy=True).transform
    geometry=transform(project,shape(feature['geometry']))
    spacing=3.5;cap=10.;gradient=.05;shore_gradient=.01
    runs=[0.,5.,10.,20.] if args.short_transitions else [0.,25.,50.,100.]
    smoothings=[25.,50.] if args.short_transitions else [0.,14.,35.]
    r,c,mask=hydrology.lake_window(geometry,terrain.grid(location,14336),4096)
    pad=int(np.ceil(6*max(smoothings)/spacing));padded=np.pad(mask,pad)
    raw=np.maximum(distance_transform_edt(padded,sampling=spacing)-spacing/2,0).astype('float32')
    distances={sigma:(np.minimum(raw,gaussian_filter(raw,sigma/spacing,mode='constant')) if sigma else raw)[pad:-pad,pad:-pad] for sigma in set(smoothings+[14.])}
    row=int(np.unravel_index(np.argmax(np.where(mask,distances[14],-1)),mask.shape)[0])
    h,w=mask.shape;x=(np.arange(w)+.5)*spacing
    def bed(distance,run):
        return cap*np.tanh(distance*gradient/cap) if run==0 else profile(distance,run,cap,shore_gradient,gradient)
    reference=bed(distances[14],0);depths={};records=[]
    for sigma in smoothings:
        for run in runs:
            z=bed(distances[sigma],run);depths[sigma,run]=z
            assert np.all((z>=0)&(z<=cap))
            assert np.all(z<=bed(distances[sigma],0)+1e-6)
            records.append(dict(shore_transition_run_m=run,distance_smoothing_sigma_m=sigma,maximum_depth_m=float(z[mask].max()),median_depth_m=float(np.median(z[mask])),area_under_1m_depth_percent=float(100*np.mean(z[mask]<1))))
    # Controlled one-bank profiles separate the designed shore shape from polygon geometry.
    d=np.linspace(0,400,4001)
    fig,axes=plt.subplots(1,2,figsize=(13,5.5),layout='constrained')
    colours=['#777777','#15918b','#7655ac','#cb7039']
    for run,colour in zip(runs,colours):
        z=bed(d,run)
        assert abs(z[0])<1e-12 and np.all(np.diff(z)>=0)
        assert np.max(np.diff(z)/np.diff(d))<=.050001
        for ax in axes:ax.plot(d,z,color=colour,lw=2,label='No shoreline transition' if run==0 else f'{run:g} m shore transition')
    axes[0].set(xlim=(0,60 if args.short_transitions else 120),ylim=(3.2 if args.short_transitions else 5.5,0),title='Shoreline close-up (profile before spatial smoothing)')
    axes[1].set(xlim=(0,400),ylim=(10.5,0),title='Shoreline to deep water')
    for ax in axes:
        ax.set(xlabel='Distance inward from bank (m)',ylabel='Depth below water surface (m)');ax.grid(alpha=.2);ax.legend(fontsize=9)
    fig.suptitle('10 m cap / 5% main gradient / gentler 1% shoreline ramping smoothly to 5%\nRun = horizontal transition distance, not a flat shelf; tanh softens slopes near the depth cap',fontsize=13)
    fig.savefig(out/'01-shoreline-profiles.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(2,len(smoothings),figsize=(6*len(smoothings),10),layout='constrained')
    shore_x=x[np.flatnonzero(mask[row])[0]]-spacing/2
    for ax,sigma in zip(axes.flat,smoothings*2):
        for run,colour in zip(runs,colours):
            ax.plot(x,np.where(mask[row],depths[sigma,run][row],np.nan),color=colour,lw=2,label=f'{run:g} m shore transition')
        ax.plot(x,np.where(mask[row],reference[row],np.nan),color='black',ls='--',lw=1,label='Previous trial: 14 m smooth / no transition')
        ax.set(title=f'{sigma:g} m distance smoothing',xlabel='Across the same lagoon section (m)',ylabel='Depth below water surface (m)',xlim=(0,w*spacing),ylim=(10.5,0));ax.grid(alpha=.2);ax.legend(fontsize=8)
    for ax in axes[1]:
        ax.set(xlim=(shore_x,shore_x+60),ylim=(3.2,0),title=ax.get_title()+' / shore close-up')
    fig.suptitle('Actual Waiuna cross-section / top row: 0–10.5 m depth; bottom row: 0–3.2 m depth\nCap 10 m · main gradient 5% · initial shore gradient 1% for nonzero runs',fontsize=14)
    fig.savefig(out/'02-lagoon-cross-sections.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(len(smoothings),len(runs),figsize=(15,5*len(smoothings)),layout='constrained')
    for i,sigma in enumerate(smoothings):
        for j,run in enumerate(runs):
            z=depths[sigma,run];rec=records[i*len(runs)+j];ax=axes[i,j]
            im=ax.imshow(np.where(mask,z,np.nan),cmap='Blues',vmin=0,vmax=10,extent=(0,w*spacing,h*spacing,0))
            ax.axhline((row+.5)*spacing,color='#c56a18',ls='--',lw=.6)
            ax.set(title=f'{run:g} m shore run / {sigma:g} m smoothing\nMax {rec["maximum_depth_m"]:.1f} m · median {rec["median_depth_m"]:.1f} m',xlabel='Map x (m)',ylabel='Map y (m)')
    fig.colorbar(im,ax=axes,label='Depth below water surface (m) — shared scale',shrink=.6)
    fig.suptitle('Waiuna / columns vary shoreline run; rows vary distance smoothing\nSame polygon, 10 m cap and 5% main gradient throughout · dashed line: common section',fontsize=14)
    fig.savefig(out/'03-depth-maps.png',dpi=160);plt.close(fig)
    result=dict(cap_m=cap,main_gradient=gradient,initial_shore_gradient=shore_gradient,pixel_spacing_m=spacing,ideal_bank_distance_to_1m_depth_m={str(run):float(np.interp(1,bed(d,run),d)) for run in runs},source_polygon=str(source.relative_to(root)),shore_profile='Integrated smoothstep gradient transition, followed by tanh depth cap',smoothing='Existing lagoon distance smoothing: min(raw distance, Gaussian-filtered distance); sigma in metres. Not a shallow shelf or river crease diffusion.',results=records)
    (out/'parameters-and-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(records,indent=2))


if __name__=='__main__':main()
