"""Nine candidates: 1/2/3m bank transitions × 1/2/3m local smoothing.
Run with uv run --all-packages python packages/water-carver/examples/plot_finalists.py
"""
import csv
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from plot_two_pass import profile, smooth_creases
from plot_samples import ROOT

OUT=ROOT/'output/water-carver-samples/two-pass-finalists'
WIDTHS=[5,10,20,40,80,160,400]
TRANSITIONS=[1,2,3]
SMOOTHING=[1,2,3]
COLORS=['#007f86','#7a52a4','#c17a20']
OUT.mkdir(parents=True,exist_ok=True)
curves={};records=[]
for width in WIDTHS:
    x=np.linspace(0,width,int(width/.1)+1);d=np.minimum(x,width-x)
    for transition in TRANSITIONS:
        raw=profile(d,transition)
        for sigma in SMOOTHING:
            z,_=smooth_creases(raw,d,.1,sigma)
            assert np.isfinite(z).all() and z.min()>=-1e-10 and z.max()<=10+1e-10
            np.testing.assert_array_equal(z[d<=.5],raw[d<=.5])
            np.testing.assert_allclose(z,z[::-1],atol=1e-10)
            assert np.all(np.diff(z[:len(z)//2+1])>=-1e-10)
            curves[width,transition,sigma]=(x,raw,z)
            records.append(dict(width_m=width,transition_m=transition,smoothing_m=sigma,raw_max_depth_m=float(raw.max()),max_depth_m=float(z.max()),depth_loss_m=float(raw.max()-z.max())))

def width_plot(widths,path):
    fig,axes=plt.subplots(len(widths),3,figsize=(17,3.15*len(widths)+1),layout='constrained',squeeze=False)
    for row,width in enumerate(widths):
        limit=max(np.max(curves[width,t,1][1]) for t in TRANSITIONS)*1.1
        for col,t in enumerate(TRANSITIONS):
            ax=axes[row,col]
            x,raw,_=curves[width,t,1]
            ax.plot(x,-raw,color='#8d9298',ls='--',lw=1.5,label='Before smoothing')
            for sigma,color in zip(SMOOTHING,COLORS):
                _,_,z=curves[width,t,sigma]
                ax.plot(x,-z,color=color,lw=2,label=f'{sigma} m smoothing')
            ax.set(title=f'{width} m channel · {t} m bank transition',xlabel='Across channel (m)',ylabel='Bed relative to water (m)',ylim=(-limit,limit*.035))
            ax.grid(alpha=.2)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=4)
    fig.suptitle('Nine finalists / 1–3 m bank transitions × 1–3 m smoothing\nBank-distance method · depth ceiling 10 m · same vertical scale within each row',fontsize=17)
    fig.savefig(OUT/path,dpi=170);plt.close(fig)

width_plot(WIDTHS,'01-all-seven-widths.png')
width_plot([5,10,20],'02-small-channels.png')
# Annotated numeric matrix: candid about loss of useful depth.
fig,axes=plt.subplots(1,3,figsize=(13,5),layout='constrained')
for ax,width in zip(axes,[5,10,20]):
    values=np.array([[curves[width,t,s][2].max() for s in SMOOTHING] for t in TRANSITIONS])
    im=ax.imshow(values,cmap='YlGnBu',vmin=0,vmax=4.5)
    for r in range(3):
        for c in range(3):ax.text(c,r,f'{values[r,c]:.2f} m',ha='center',va='center',fontsize=13,color='white' if values[r,c]>2.5 else '#163343')
    ax.set(title=f'{width} m-wide channel',xticks=range(3),xticklabels=['1 m','2 m','3 m'],yticks=range(3),yticklabels=['1 m','2 m','3 m'],xlabel='Smoothing scale',ylabel='Bank transition')
fig.colorbar(im,ax=list(axes),shrink=.7,label='Maximum generated depth (m)')
fig.suptitle('Small-channel depth / measured after both passes',fontsize=18)
fig.savefig(OUT/'03-depth-matrix.png',dpi=170);plt.close(fig)
with (OUT/'depth-comparison.csv').open('w',newline='') as handle:
    writer=csv.DictWriter(handle,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
(OUT/'parameters.json').write_text(json.dumps(dict(widths_m=WIDTHS,bank_transitions_m=TRANSITIONS,smoothing_scales_m=SMOOTHING,depth_cap_m=10,shore_gradient=.05,interior_gradient=.5,shore_lock_m=.5,shore_blend_m=1.5,sample_spacing_m=.1,records=records),indent=2)+'\n')
(OUT/'README.md').write_text('''# Nine two-pass finalists

Bank transitions of **1, 2, 3 m**, crossed with smoothing scales of **1, 2, 3 m**. “Bank transition” is horizontal distance inward from the bank, not target water depth. Every case retains the same 10 m depth ceiling, 5% initial shore slope, 50% nominal deeper-bank slope, 0.5 m shore lock and 1.5 m smoothing fade near shore. No width-dependent parameters, production changes or map edits.

## Small channels — start here

![Small channels](02-small-channels.png)

## All seven widths

5, 10, 20, 40, 80, 160 and 400 m. Depth axes match across each row but differ between rows. Grey dashed curves are before smoothing; teal/purple/orange are 1/2/3 m smoothing.

![All widths](01-all-seven-widths.png)

## Resulting depths

![Depth matrix](03-depth-matrix.png)

Exact results for all 63 cases: [CSV](depth-comparison.csv).

These are ideal cross-sections sampled at 0.1 m for shape comparison, not game-raster exports. Features smaller than the existing 3.5 m terrain grid will need export-resolution review. Smoothing rounds deepest points by filling them; it does not preserve depth or guarantee water flow. The broad-channel plateau is unchanged by this parameter selection.

Reproduce from repo root:

```sh
uv run --all-packages python packages/water-carver/examples/plot_finalists.py
```
''')
for width in [5,10,20]:
    print(width,[[round(float(curves[width,t,s][2].max()),3) for s in SMOOTHING] for t in TRANSITIONS])
