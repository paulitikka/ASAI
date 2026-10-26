"""Plot the four Final12 figures using the manuscript tables as the baseline.

Run from an extracted ASAI release: python plot_final12_figures.py
This script does not recompute or overwrite any manuscript table.
"""
from pathlib import Path
import csv
import json
from io import BytesIO
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parent
BLUE, ORANGE, GREY = '#0072B2', '#D55E00', '#777777'
AA = list('ACDEFGHIKLMNPQRSTVWY')

def read_csv(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))

def number(value):
    return float(value.replace('\u2212', '-').replace('%', '').strip())

def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,
                         'axes.spines.top':False,'axes.spines.right':False,
                         'axes.linewidth':0.8,'pdf.fonttype':42,'svg.fonttype':'none'})
    out = ROOT/'manuscript_figures'
    out.mkdir(exist_ok=True)
    def save(fig, stem):
        for ext in ['pdf','svg','png']:
            buffer=BytesIO()
            fig.savefig(buffer, format=ext, dpi=600, facecolor='white')
            destination=out/(stem+'.'+ext)
            temporary=destination.with_suffix('.'+ext+'.tmp')
            temporary.write_bytes(buffer.getvalue())
            os.replace(temporary,destination)
        plt.close(fig)

    # Figure 1 uses the already verified residue histogram counts.
    hist = read_csv(ROOT/'source_data/Figure1_histogram_bins.csv')
    fig, axes = plt.subplots(4,5,figsize=(7.2,7.8),layout='constrained')
    for ax, aa in zip(axes.flat, AA):
        rows=[r for r in hist if r['aa']==aa]
        ax.bar([float(r['bin_left']) for r in rows],
               [float(r['density']) for r in rows],width=0.025,align='edge',
               color=BLUE,linewidth=0)
        ax.set_title(aa,loc='left',weight='bold',fontsize=10)
        ax.set_xlim(0,1.5)
        ax.set_xticks([0,.5,1,1.5])
        ax.tick_params(labelsize=7)
    fig.supxlabel('Residue relative solvent accessibility')
    fig.supylabel('Probability density')
    save(fig,'Figure_1_residue_RSA')

    # Figure 2 combines frozen-data SDs with the PC2 coordinates printed in Table 1.
    table1={r['Type']:r for r in read_csv(ROOT/'manuscript_tables/Table_1.csv')}
    sd={r['aa']:r for r in read_csv(ROOT/'source_data/Figure2_PC2_dispersion.csv')}
    x=np.array([float(sd[a]['sd_chain_mean_RSA']) for a in AA])
    y=np.array([number(table1[a]['PC2']) for a in AA])
    fit=np.polyfit(x,y,1)
    r,p=pearsonr(x,y); rho,sp=spearmanr(x,y)
    fig, ax=plt.subplots(figsize=(5.5,4.1),layout='constrained')
    ax.scatter(x,y,color='#252525',s=17,zorder=3)
    xx=np.array([x.min()-.001,x.max()+.001])
    ax.plot(xx,np.polyval(fit,xx),color=ORANGE,lw=1.5)
    dx=np.array([0,-.038044,.042,.040338,-.040338,.042,-.036418,-.021389,.034408,.005852,
                 -.040338,.036418,.034408,0,-.0167,-.021389,-.042,0,-.042,-.031955])
    dy=np.array([-.06,-.025420,0,.016709,.016709,0,.029888,-.051637,.034408,-.059415,
                 .016709,-.029888,.034408,.06,.055053,-.051637,0,-.06,0,-.038937])
    normal=np.array([-fit[0]*np.ptp(x)/np.ptp(y),1.])
    normal/=np.linalg.norm(normal)
    offsets=np.column_stack([dx,dy])
    for i in range(20):
        side=np.sign(y[i]-np.polyval(fit,x[i]))
        component=offsets[i]@normal
        if side*component<0:offsets[i]-=2*component*normal
    for i,(a, xi, yi) in enumerate(zip(AA,x,y)):
        lab=(xi+offsets[i,0]*np.ptp(x),yi+offsets[i,1]*np.ptp(y))
        lab={'A':(.0945,-.010),'T':(.0905,.050),
             'S':(.1010,-.070),'Y':(.0965,-.046)}.get(a,lab)
        ax.annotate(a,(xi,yi),xytext=lab,
                    ha='center',va='center',fontsize=8,color='white',weight='bold',
                    bbox={'boxstyle':'circle,pad=.2','fc':BLUE,'ec':'#252525','lw':.7},
                    arrowprops={'arrowstyle':'-','color':GREY,'lw':.6},zorder=4)
    ax.text(.98,.98,f'Pearson r = {r:.3f}',transform=ax.transAxes,
            ha='right',va='top',weight='bold')
    ax.set_xlabel('SD of per-chain mean RSA')
    ax.set_ylabel('PC2 coordinate')
    ax.set_xlim(.062,.158); ax.set_ylim(-.67,.43)
    save(fig,'Figure_2_PC2_dispersion')
    with (out/'Figure_2_statistics.json').open('w') as handle:
        json.dump({'n_types':20,'Pearson_r':float(r),'Pearson_P':float(p),
                   'Spearman_rho':float(rho),'Spearman_P':float(sp),
                   'PC2_source':'manuscript_tables/Table_1.csv',
                   'SD_source':'source_data/Figure2_PC2_dispersion.csv',
                   'interpretation':'Descriptive; both quantities use the same chain distributions.'},handle,indent=2)

    # Figure 3 uses the printed CV values without substituting a new analysis run.
    s1=read_csv(ROOT/'manuscript_tables/Table_S1.csv')
    cv=[100*number(r['CV across types']) for r in s1]
    labels=['0']+[f'{i/10:g}\u2013{(i+1)/10:g}' for i in range(9)]+['>0.9']
    fig,ax=plt.subplots(figsize=(7.2,4.3),layout='constrained')
    ax.plot(range(11),cv,color=BLUE,lw=1.8)
    ax.scatter(range(11),cv,s=42,color=ORANGE,edgecolor='#252525',zorder=3)
    ax.set_xticks(range(11),labels,fontsize=8)
    ax.set_xlabel('Per-chain mean RSA interval')
    ax.set_ylabel('Coefficient of variation (%)')
    ax.set_ylim(0,max(cv)*1.1)
    save(fig,'Figure_3_divergence')

    # Figure 4 preserves the paired-point design of the user-selected PDF.
    s2={r['Type']:r for r in read_csv(ROOT/'manuscript_tables/Table_S2.csv')}
    fig,axes=plt.subplots(1,2,figsize=(7.2,5.4),sharey=True,layout='constrained')
    aa_order=AA[::-1]
    for j,(ax,pc,share) in enumerate(zip(axes,['PC1','PC2'],[96.48,3.23])):
        e=np.array([number(s2[a][pc+' euk']) for a in aa_order])
        p_=np.array([number(s2[a][pc+' prok']) for a in aa_order])
        yy=np.arange(20)
        ax.hlines(yy,p_,e,color=GREY,lw=.9)
        ax.scatter(e,yy,s=25,facecolor=BLUE,edgecolor='#252525',lw=.7,label='Eukaryotes',zorder=3)
        ax.scatter(p_,yy,s=25,facecolor='white',edgecolor=ORANGE,lw=1.1,label='Prokaryotes',zorder=3)
        ax.axvline(0,color='#aaaaaa',ls=':',lw=.8,zorder=0)
        ax.set_yticks(yy,aa_order)
        ax.set_ylim(-.7,19.7)
        ax.set_xlabel(f'{pc} coordinate\nVariance explained: {share:.2f}%')
        ax.set_title('ab'[j],loc='left',weight='bold')
    axes[0].set_ylabel('Amino-acid type')
    axes[1].legend(loc='lower right',bbox_to_anchor=(1,1.02),ncol=2,frameon=False,fontsize=8)
    save(fig,'Figure_4_kingdom_coordinates')

if __name__=='__main__':
    main()
