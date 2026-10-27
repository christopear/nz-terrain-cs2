"""Steeper-bank trial: 100% nominal bank gradient; all plots share a depth scale.
Run with uv run --all-packages python packages/water-carver/examples/plot_steeper_banks.py
"""
import csv
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from plot_two_pass import profile, smooth_creases
from plot_samples import ROOT

OUT=ROOT/'output/water-carver-samples/two-pass-steeper-banks'
WIDTHS=[5,10,20,40,80,160,400]
TRANSITIONS=[1,2,3]
SMOOTHING=[1,2,3]
COLORS=['#007f86','#7a52a4','#c17a20']
OUT.mkdir(parents=True,exist_ok=True)
curves={};records=[]
for width in WIDTHS:
    x=np.linspace(0,width,int(width/.1)+1);d=np.minimum(x,width-x)
    for transition in TRANSITIONS:
        raw=profile(d,transition,deep_gradient=1.0)
        for sigma in SMOOTHING:
            z,_=smooth_creases(raw,d,.1,sigma)
            assert np.isfinite(z).all() and z.min()>=-1e-10 and z.max()<=10+1e-10
            np.testing.assert_array_equal(z[d<=.5],raw[d<=.5])
            np.testing.assert_allclose(z,z[::-1],atol=1e-10)
            assert np.all(np.diff(z[:len(z)//2+1])>=-1e-10)
            curves[width,transition,sigma]=(x,raw,z)
            records.append(dict(width_m=width,transition_m=transition,smoothing_m=sigma,raw_max_depth_m=float(raw.max()),max_depth_m=float(z.max()),depth_loss_m=float(raw.max()-z.max()),maximum_slope_percent=float(np.max(np.abs(np.gradient(z,.1)))*100)))

def width_plot(widths,path):
    fig,axes=plt.subplots(len(widths),3,figsize=(17,3.15*len(widths)+1),layout='constrained',squeeze=False)
    for row,width in enumerate(widths):
        limit=10.5
        for col,t in enumerate(TRANSITIONS):
            ax=axes[row,col]
            x,raw,_=curves[width,t,1]
            ax.plot(x,-raw,color='#8d9298',ls='--',lw=1.5,label='Before smoothing')
            for sigma,color in zip(SMOOTHING,COLORS):
                _,_,z=curves[width,t,sigma]
                ax.plot(x,-z,color=color,lw=2,label=f'{sigma} m smoothing')
            depth_labels=' / '.join(f'{curves[width,t,s][2].max():.2f}' for s in SMOOTHING)
            ax.set(title=f'{width} m channel · {t} m bank transition\nDepths: {depth_labels} m (1 / 2 / 3 m smoothing)',xlabel='Across channel (m)',ylabel='Bed relative to water (m)',ylim=(-10.5,.5),yticks=[0,-2,-4,-6,-8,-10])
            ax.grid(alpha=.2)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=4)
    fig.suptitle('Steeper banks / 1–3 m transitions × 1–3 m smoothing\n100% nominal bank slope · 10 m ceiling · IDENTICAL depth scale on every panel',fontsize=17)
    fig.savefig(OUT/path,dpi=170);plt.close(fig)

width_plot(WIDTHS,'01-all-seven-widths.png')
width_plot([10,20,40],'02-key-channel-widths.png')
# Annotated numeric matrix: candid about loss of useful depth.
fig,axes=plt.subplots(1,4,figsize=(16,5),layout='constrained')
for ax,width in zip(axes,[5,10,20,40]):
    values=np.array([[curves[width,t,s][2].max() for s in SMOOTHING] for t in TRANSITIONS])
    im=ax.imshow(values,cmap='YlGnBu',vmin=0,vmax=10)
    for r in range(3):
        for c in range(3):ax.text(c,r,f'{values[r,c]:.2f} m',ha='center',va='center',fontsize=13,color='white' if values[r,c]>6 else '#163343')
    ax.set(title=f'{width} m-wide channel',xticks=range(3),xticklabels=['1 m','2 m','3 m'],yticks=range(3),yticklabels=['1 m','2 m','3 m'],xlabel='Smoothing scale',ylabel='Bank transition')
fig.colorbar(im,ax=list(axes),shrink=.7,label='Maximum generated depth (m)')
fig.suptitle('Achieved depths / measured after both passes',fontsize=18)
fig.savefig(OUT/'03-depth-matrix.png',dpi=170);plt.close(fig)
with (OUT/'depth-comparison.csv').open('w',newline='') as handle:
    writer=csv.DictWriter(handle,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
(OUT/'parameters.json').write_text(json.dumps(dict(widths_m=WIDTHS,bank_transitions_m=TRANSITIONS,smoothing_scales_m=SMOOTHING,depth_cap_m=10,shore_gradient=.05,interior_gradient=1.0,shore_lock_m=.5,shore_blend_m=1.5,sample_spacing_m=.1,records=records),indent=2)+'\n')
# Focused target check: keep 2m transition + 1m smoothing; compare depth ceilings.
cap_records=[]
fig,axes=plt.subplots(2,4,figsize=(16,9),layout='constrained')
for ax,width in zip(axes.flat,WIDTHS):
    x=np.linspace(0,width,int(width/.1)+1);d=np.minimum(x,width-x)
    for cap,color in [(8,'#007f86'),(10,'#7a52a4')]:
        raw=profile(d,2,cap=cap,deep_gradient=1.)
        z,_=smooth_creases(raw,d,.1,1)
        assert np.isfinite(z).all() and z.max()<=cap+1e-10
        assert np.all(np.diff(z[:len(z)//2+1])>=-1e-10)
        cap_records.append(dict(width_m=width,cap_m=cap,transition_m=2,smoothing_m=1,max_depth_m=float(z.max())))
        ax.plot(x,-z,color=color,lw=2,label=f'{cap} m cap → {z.max():.2f} m deep')
    ax.set(title=f'{width} m-wide channel',xlabel='Across channel (m)',ylabel='Bed elevation (m)',ylim=(-10.5,.5),yticks=[0,-2,-4,-6,-8,-10]);ax.legend(fontsize=9);ax.grid(alpha=.2)
axes.flat[-1].axis('off')
axes.flat[-1].text(.05,.85,'Fixed in this comparison:\n\n2 m shoreline transition\n1 m smoothing scale\n100% nominal bank gradient\n\nOnly the depth ceiling changes.\nAll panels have the same depth scale.',va='top',fontsize=12)
fig.suptitle('Target check / 8 m versus 10 m depth ceiling',fontsize=18)
fig.savefig(OUT/'04-eight-versus-ten-metre-cap.png',dpi=170);plt.close(fig)
(OUT/'cap-comparison.json').write_text(json.dumps(cap_records,indent=2)+'\n')
print('CAP COMPARISON',json.dumps(cap_records))
(OUT/'README.md').write_text("""# Steeper banks: nine candidates

This keeps the bank-distance approach and local crease smoothing. The only modelling change from the previous nine-candidate trial is the nominal deeper-bank gradient: **100% (45°)** instead of 50% (~26.6°). The initial shore gradient remains 5%; tanh capping reduces the realised bank slope below the nominal value. No width normalisation or mandatory stream depth is introduced.

All plots use **the same −10.5 to +0.5 m vertical scale**. Horizontal scales vary with channel width, so apparent line angles are not physical bank angles. Every panel labels the maximum depths after 1 / 2 / 3 m smoothing.

## Target check: 8 m versus 10 m ceiling

Fixed 2 m bank transition and 1 m smoothing; only the depth ceiling changes.

![Cap comparison](04-eight-versus-ten-metre-cap.png)

## Key widths: 10, 20, 40 m

![Key widths](02-key-channel-widths.png)

## All requested widths

![All widths](01-all-seven-widths.png)

## Numeric comparison

![Depth matrix](03-depth-matrix.png)

Exact depths, depth loss and measured maximum slopes for all 63 cases: [CSV](depth-comparison.csv).

Parameters: 1/2/3 m bank transitions × 1/2/3 m smoothing; 10 m ceiling; initial slope 5%; nominal interior slope 100%; shore lock 0.5 m; shore blend 1.5 m; 0.1 m cross-section sampling. The smoothing remains the previous filling/averaging method, not an excavation-only method.

This is a shape trial, not an in-game flow result. A 3.5 m export raster will not retain all details of these 1–3 m transitions. No production defaults or maps changed.

Reproduce:

```sh
uv run --all-packages python packages/water-carver/examples/plot_steeper_banks.py
```
""")
for width in WIDTHS:
    print(width,[[round(float(curves[width,t,s][2].max()),3) for s in SMOOTHING] for t in TRANSITIONS])
