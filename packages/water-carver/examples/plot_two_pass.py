"""Bank-driven S-profile and local crease smoothing experiments.

Run: uv run --all-packages python packages/water-carver/examples/plot_two_pass.py
No local/global channel widths are input to either algorithm. Widths below only
construct test polygons/sections. Core package and map outputs are unchanged.
"""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import distance_transform_edt, map_coordinates
from shapely.geometry import LineString
from shapely.ops import unary_union
from water_carver import BedParameters, generate_bed
from plot_samples import ROOT, grid_for

OUT=ROOT/'output/water-carver-samples/two-pass-banks'
WIDTHS=[5,10,20,40,80,160,400]
VARIANTS=[('5 m bank transition',5),('10 m bank transition',10),('Compact 2 m transition',2)]


from water_carver.river import profile, smooth_creases

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    # One-dimensional channel sections: same routines for every width, metre units.
    records=[]
    fig,axes=plt.subplots(7,3,figsize=(18,23),layout='constrained')
    for row,width in enumerate(WIDTHS):
        x=np.linspace(0,width,int(width/.1)+1);spacing=x[1]-x[0]
        d=np.minimum(x,width-x)
        rawmax=max(np.max(profile(d,t)) for _,t in VARIANTS)
        for col,(label,transition) in enumerate(VARIANTS):
            raw=profile(d,transition);ax=axes[row,col]
            ax.plot(x,-raw,color='#939393',ls='--',lw=1.7,label='Pass 1: banks only')
            depths={'raw_max_m':float(raw.max())}
            for sigma,color in [(2,'#007f86'),(5,'#b66c2c')]:
                smoothed,weight=smooth_creases(raw,d,spacing,sigma)
                assert np.isfinite(smoothed).all() and smoothed.min()>=-1e-10 and smoothed.max()<=10+1e-10
                np.testing.assert_array_equal(smoothed[d<=.5],raw[d<=.5])
                np.testing.assert_allclose(smoothed,smoothed[::-1],atol=1e-10)
                assert np.all(np.diff(smoothed[:len(smoothed)//2+1])>=-1e-10), "Unexpected central hump"
                ax.plot(x,-smoothed,color=color,lw=2,label=f'Pass 2: {sigma} m smoothing scale')
                depths[f'smoothed_{sigma}m_max_m']=float(smoothed.max())
            records.append(dict(width_m=width,bank_transition_m=transition,**depths))
            ax.set(title=f'{width} m channel · {label}',xlabel='Across channel (m)',ylabel='Bed elevation (m)',ylim=(-max(.5,rawmax*1.12),.03*max(.5,rawmax)))
            ax.grid(alpha=.18)
            if row==0:ax.legend(fontsize=8,loc='lower center')
    fig.suptitle('Two passes / bank-shaped depth, then local bottom rounding\nNo channel-width input · 10 m depth ceiling · rows use different depth scales to show small-stream detail',fontsize=18)
    fig.savefig(OUT/'01-all-widths.png',dpi=160);plt.close(fig)

    # Bank-side controls and numerical depth trade-off, kept readable separately.
    fig,axes=plt.subplots(1,2,figsize=(14,6),layout='constrained')
    d=np.linspace(0,80,801)
    for (label,t),color in zip(VARIANTS,['#007f86','#7a52a4','#b66c2c']):
        axes[0].plot(d,-profile(d,t),color=color,lw=2,label=label)
    axes[0].set(title='Pass 1: gentle edge → steeper bank → capped depth',xlabel='Distance from nearest bank (m)',ylabel='Bed relative to water (m)');axes[0].legend();axes[0].grid(alpha=.2)
    rows=[r for r in records if r['bank_transition_m']==5]
    for key,label,color in [('raw_max_m','Banks only','#939393'),('smoothed_2m_max_m','2 m smoother','#007f86'),('smoothed_5m_max_m','5 m smoother','#b66c2c')]:
        axes[1].plot(WIDTHS,[r[key] for r in rows],marker='o',color=color,label=label)
    axes[1].set(title='Deepest point · 5 m bank transition',xlabel='Channel width (m; logarithmic scale)',ylabel='Maximum depth (m)',xscale='log',xticks=WIDTHS,xticklabels=WIDTHS);axes[1].legend();axes[1].grid(alpha=.2)
    fig.suptitle('The depth trade-off is explicit: smoothing rounds the floor by making it shallower',fontsize=16)
    fig.savefig(OUT/'02-bank-shapes-and-depths.png',dpi=170);plt.close(fig)

    # 2D test: main channel and tributary, island, no global width normalisation.
    main=LineString([(0,90),(80,90),(160,110),(240,90),(320,90)]).buffer(20,cap_style='flat',quad_segs=16)
    branch=LineString([(155,110),(175,170),(220,220)]).buffer(8,cap_style='flat',quad_segs=16)
    island=LineString([(72,88),(86,88)]).buffer(4,quad_segs=16)
    water=unary_union([main,branch]).difference(island)
    grid=grid_for((-5,55,325,235),1)
    base=generate_bed(water,grid,water_level_m=0,parameters=BedParameters(1000,1,'linear'))
    raw=profile(base.depth_m,5)
    smoothed,weight=smooth_creases(raw,base.depth_m,1,5)
    np.testing.assert_array_equal(np.isfinite(raw),np.isfinite(smoothed))
    assert np.nanmin(smoothed)>=0 and np.nanmax(smoothed)<=10
    # A fully disconnected channel cannot influence a neighbour through dry cells.
    a=np.array([0.,2.,3.,np.nan,0.,8.,9.]);dist=np.array([0.,1.,2.,np.nan,0.,1.,2.])
    b=a.copy();b[5:]=100
    sa,_=smooth_creases(a,dist,1,2);sb,_=smooth_creases(b,dist,1,2)
    np.testing.assert_array_equal(sa[:3],sb[:3])
    fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained')
    aff=grid.transform;h,w=grid.shape;extent=(aff.c,aff.c+w,aff.f-h,aff.f)
    for ax,values,title,cm,lo,hi in [(axes[0,0],raw,'Pass 1: banks only','YlGnBu',0,10),(axes[0,1],smoothed,'Pass 2: 5 m local smoothing','YlGnBu',0,10),(axes[1,0],smoothed-raw,'Depth change (negative = shallower)','RdBu',-1.5,1.5),(axes[1,1],np.where(base.mask,weight,np.nan),'Local diffusion strength','magma',0,1)]:
        if 'change' in title:
            bound=max(.01,float(np.nanmax(np.abs(smoothed-raw))));lo=-bound;hi=bound
        im=ax.imshow(values,extent=extent,origin='upper',vmin=lo,vmax=hi,cmap=cm,interpolation='nearest')
        ax.set(title=title,xlabel='Easting (m)',ylabel='Northing (m)');fig.colorbar(im,ax=ax,shrink=.7,label='Blend strength' if 'diffusion strength' in title else 'Metres')
    fig.suptitle('2D junction test / main channel + tributary + island\nBank distances only; smoothing exchanges depth between adjacent water cells, never across land',fontsize=17)
    fig.savefig(OUT/'03-junction-test.png',dpi=170);plt.close(fig)
    (OUT/'parameters.json').write_text(json.dumps(dict(depth_cap_m=10,shore_gradient=.05,interior_gradient=.5,bank_transitions_m=[5,10,2],smoothing_scales_m=[2,5],shore_lock_m=.5,shore_blend_m=1.5,crease_influence_radius='2 × smoothing scale',depths=records,junction_max_absolute_change_m=float(np.nanmax(np.abs(smoothed-raw)))),indent=2)+'\n')
    (OUT/'README.md').write_text('''# Two-pass bank-driven experiments

    No local or global channel width is supplied to either pass. The example widths only create test channels. No production carver behaviour or maps have changed.

    ## All requested widths

    5, 10, 20, 40, 80, 160 and 400 m. Each row uses a common depth scale across its three columns, but scales differ between rows so tiny-stream behaviour stays visible. Grey = unsmoothed; teal = 2 m smoother; orange = 5 m smoother.

    ![All widths](01-all-widths.png)

    ## Bank profiles and the depth cost

    ![Profiles and depths](02-bank-shapes-and-depths.png)

    ## Two-dimensional tributary/island test

    ![Junction test](03-junction-test.png)

    ### Pass 1

    Start at a 5% shoreline slope, increase smoothly to a nominal 50% bank slope over 5 or 10 m (plus a compact 2 m comparison), then use tanh to approach a 10 m ceiling. The integrated slope transition gives an S-like depth curve. Capping reduces the realised slopes below their uncapped values.

    ### Pass 2

    Identify sharp ridges in the distance-to-bank field, then apply masked diffusion only within twice the selected smoothing scale of those ridges. Diffusion strength fades toward the edge of this region; it is not a post-processing blend. Cells within 0.5 m of shore are fixed; the next 1.5 m gradually permits changes. Smoothing passes values only between adjacent water cells, never through land. It also rounds geometric corners and junction creases; it is not a classified river-bottom detector. The bank profile in a very narrow stream can be altered beyond the protected 0.5 m edge because its bottom is also close to shore.

    Smoothing generally fills deepest points, reducing maximum depth. There is no depth or volume preservation, no guarantee of a 10–20 m flat floor, no forced 10 m stream depth, and no flow simulation. Very narrow streams remain shallow; 5–10 m gentle margins cannot coexist with deep flat floors inside 5–10 m channels. The 400 m example still approaches a broad 10 m plateau: local crease smoothing alone cannot redesign a flat region.

    These are high-resolution mathematical examples. At the existing map's 3.5 m raster spacing, 5 m channels span only about 1–2 cells, so the fine bank shape would not be retained.

    Reproduce:

    ```sh
    uv run --all-packages python packages/water-carver/examples/plot_two_pass.py
    ```
    ''')
    print(json.dumps({'five_metre_transition_depths':rows,'junction_max_change_m':float(np.nanmax(np.abs(smoothed-raw)))},indent=2))


if __name__ == "__main__":
    main()
