"""Compare synthetic Waiuna lagoon beds; leaves map packages/settings unchanged."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from pyproj import Transformer
from scipy.ndimage import distance_transform_edt, gaussian_filter
from shapely.geometry import shape
from shapely.ops import transform
from nz_terrain import terrain, hydrology


def main():
    root=terrain.BASE
    out=root/'output/jamestown-river-lab/waiuna/variations';out.mkdir(parents=True,exist_ok=True)
    source=root/'output/jamestown-polygon-rivers/metadata/lake-polygons.geojson'
    feature=next(f for f in json.loads(source.read_text())['features'] if f['properties'].get('name')=='Waiuna Lagoon')
    location=json.loads((root/'locations.json').read_text())['jamestown']
    project=Transformer.from_crs(4326,terrain.projection(location['longitude'],location['latitude']),always_xy=True).transform
    geometry=transform(project,shape(feature['geometry']))
    spacing=3.5;smoothing=14.;level=15.
    r,c,original=hydrology.lake_window(geometry,terrain.grid(location,14336),4096)
    pad=int(np.ceil(6*smoothing/spacing));mask=np.pad(original,pad)
    distance=np.maximum(distance_transform_edt(mask,sampling=spacing)-spacing/2,0).astype('float32')
    distance=np.minimum(distance,gaussian_filter(distance,smoothing/spacing,mode='constant'))
    distance=distance[pad:-pad,pad:-pad]
    row=int(np.unravel_index(np.argmax(np.where(original,distance,-1)),original.shape)[0])
    h,w=original.shape;x=(np.arange(w)+.5)*spacing
    caps=[5.,10.,15.];gradients=[.02,.05,.10];records=[];depths={}
    for cap in caps:
        for gradient in gradients:
            depth=cap*np.tanh(distance/(cap/gradient))
            check=hydrology.basin_depth(mask,spacing,cap,cap/gradient,smoothing)[pad:-pad,pad:-pad]
            np.testing.assert_allclose(depth,check)
            assert np.all((depth>=0)&(depth<=cap))
            depths[cap,gradient]=depth
            values=depth[original]
            records.append(dict(cap_m=cap,gradient_percent=gradient*100,maximum_depth_m=float(values.max()),median_depth_m=float(np.median(values)),area_at_least_2m_deep_percent=float(100*np.mean(values>=2)),deepest_bed_elevation_m=level-float(values.max())))
    fig,axes=plt.subplots(3,3,figsize=(12,16),layout='constrained')
    for i,cap in enumerate(caps):
        for j,gradient in enumerate(gradients):
            ax=axes[i,j];z=depths[cap,gradient];rec=records[i*3+j]
            im=ax.imshow(np.where(original,z,np.nan),cmap='Blues',vmin=0,vmax=15,extent=(0,w*spacing,h*spacing,0))
            ax.axhline((row+.5)*spacing,color='#c56a18',ls='--',lw=.8)
            current=' · CURRENT' if cap==10 and gradient==.02 else ''
            ax.set(title=f'{cap:g} m cap · {gradient*100:g}% gradient{current}\nDeepest {rec["maximum_depth_m"]:.1f} m · median {rec["median_depth_m"]:.1f} m',xlabel='Along map x (m)',ylabel='Along map y (m)')
    fig.colorbar(im,ax=axes,label='Depth below lagoon water surface (m) — identical scale for every panel',shrink=.6)
    fig.suptitle('Waiuna Lagoon / same polygon, different synthetic beds\nRows: depth ceiling · Columns: shore-gradient parameter · dashed line: shared cross-section',fontsize=15)
    fig.savefig(out/'01-depth-maps.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(16,9),layout='constrained')
    colours=['#15918b','#7655ac','#cb7039']
    for i,cap in enumerate(caps):
        ax=axes[0,i]
        for gradient,colour in zip(gradients,colours):
            z=np.where(original[row],depths[cap,gradient][row],np.nan)
            ax.plot(x,z,color=colour,lw=2,label=f'{gradient*100:g}% gradient')
        ax.set_title(f'Hold cap at {cap:g} m; vary gradient')
    for j,gradient in enumerate(gradients):
        ax=axes[1,j]
        for cap,colour in zip(caps,colours):
            z=np.where(original[row],depths[cap,gradient][row],np.nan)
            ax.plot(x,z,color=colour,lw=2,label=f'{cap:g} m cap')
        ax.set_title(f'Hold gradient at {gradient*100:g}%; vary cap')
    for ax in axes.flat:
        current=np.where(original[row],depths[10.,.02][row],np.nan)
        ax.plot(x,current,color='.4',ls='--',lw=1,label='Current: 10 m / 2%')
        ax.axhline(0,color='steelblue',lw=1)
        ax.set(xlabel='Across lagoon at dashed section (m)',ylabel='Depth below water surface (m)',ylim=(16,-.5),xlim=(0,w*spacing));ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Waiuna cross-sections / every panel uses the SAME 0–16 m depth scale\n0 m means the lagoon water surface (+15 m above sea level); depth increases downwards',fontsize=15)
    fig.savefig(out/'02-cross-sections.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(14,5),layout='constrained')
    for ax,(key,title,fmt,cmap,limit) in zip(axes,[('maximum_depth_m','Deepest point (m)', '.1f','Blues',15),('median_depth_m','Median depth (m)', '.1f','Blues',15),('area_at_least_2m_deep_percent','Area at least 2 m deep (%)','.0f','YlGnBu',100)]):
        values=np.array([r[key] for r in records]).reshape(3,3)
        ax.imshow(values,vmin=0,vmax=limit,cmap=cmap)
        for i in range(3):
            for j in range(3):ax.text(j,i,format(values[i,j],fmt),ha='center',va='center',color='white' if values[i,j]>limit*.6 else '#152633',fontsize=18)
        ax.set(xticks=range(3),xticklabels=['2%','5%','10%'],yticks=range(3),yticklabels=['5 m','10 m','15 m'],xlabel='Shore-gradient parameter',ylabel='Depth cap',title=title)
    fig.suptitle('Waiuna / whole-lagoon measurements (synthetic, not surveyed)',fontsize=16)
    fig.savefig(out/'03-comparison-table.png',dpi=160);plt.close(fig)
    (out/'parameters-and-results.json').write_text(json.dumps(dict(water_surface_elevation_m=level,pixel_spacing_m=spacing,distance_smoothing_m=smoothing,method='Existing lagoon profile: cap * tanh(smoothed shore distance * gradient / cap)',source_polygon=str(source.relative_to(root)),cross_section_row_in_window=row,results=records),indent=2)+'\n')
    print(json.dumps(records,indent=2))


if __name__=='__main__':main()
