"""Regenerate Figures 1-4 and their source data from the frozen RSA input.

Run: python plot_figures.py
Run run.py first to regenerate canonical full-precision tables. Only local data are used.
"""
from pathlib import Path
from io import BytesIO
import argparse
import hashlib
import json
import platform
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd
from scipy import stats
import scipy

ROOT = Path(__file__).resolve().parent
AA = list('ACDEFGHIKLMNPQRSTVWY')
NAMES = dict(zip(AA, ['Ala','Cys','Asp','Glu','Phe','Gly','His','Ile','Lys','Leu',
                     'Met','Asn','Pro','Gln','Arg','Ser','Thr','Val','Trp','Tyr']))
BLUE = '#0072B2'
ORANGE = '#D55E00'
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8,
    'axes.labelsize':8, 'axes.titlesize':9, 'xtick.labelsize':7,
    'ytick.labelsize':7, 'axes.spines.top':False, 'axes.spines.right':False,
    'axes.linewidth':0.6, 'xtick.major.width':0.6, 'ytick.major.width':0.6,
    'pdf.fonttype':42, 'ps.fonttype':42, 'svg.fonttype':'none'})

def decompose(rows):
    centered = rows - rows.mean(axis=0)
    u, s, _ = np.linalg.svd(centered, full_matrices=False)
    return u, 100 * s*s / np.sum(s*s)

def place_labels(ax, x, y, names):
    """Place twenty labels deterministically, checking text and marker overlap."""
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    bounds = ax.get_window_extent(renderer)
    point_xy = ax.transData.transform(np.column_stack([x,y]))
    occupied = []
    label_positions = {}
    choices = [(5,5),(-5,5),(5,-5),(-5,-5),(12,0),(-12,0),(0,12),(0,-12),
               (14,10),(-14,10),(14,-10),(-14,-10),(22,0),(-22,0),(0,22),(0,-22),
               (24,16),(-24,16),(24,-16),(-24,-16)]
    for i, name in enumerate(names):
        best = None
        for dx,dy in choices:
            ann = ax.annotate(name,(x[i],y[i]),xytext=(dx,dy),textcoords='offset points',
                ha='left' if dx>0 else 'right' if dx<0 else 'center',
                va='bottom' if dy>0 else 'top' if dy<0 else 'center',fontsize=7)
            bb = ann.get_window_extent(renderer).expanded(1.12,1.18)
            overlaps = sum(bb.overlaps(other) for other in occupied)
            on_points = sum(bb.contains(px,py) for px,py in point_xy)
            outside = not (bounds.contains(bb.x0,bb.y0) and bounds.contains(bb.x1,bb.y1))
            score = overlaps*10000 + on_points*1000 + outside*100000 + dx*dx+dy*dy
            ann.remove()
            if best is None or score<best[0]: best=(score,dx,dy)
        _,dx,dy = best
        ann=ax.annotate(name,(x[i],y[i]),xytext=(dx,dy),textcoords='offset points',
            ha='left' if dx>0 else 'right' if dx<0 else 'center',
            va='bottom' if dy>0 else 'top' if dy<0 else 'center',fontsize=7,
            arrowprops={'arrowstyle':'-','lw':0.35,'color':'#888888'}
                if abs(dx)+abs(dy)>20 else None)
        bb=ann.get_window_extent(renderer).expanded(1.08,1.12)
        occupied.append(bb)
        # Text centre coordinates also allow R to reproduce label placement.
        lx,ly=ax.transData.inverted().transform([(bb.x0+bb.x1)/2,(bb.y0+bb.y1)/2])
        label_positions[name]=(float(lx),float(ly))
    return label_positions

def save(fig, stem, figures, combined):
    for extension in ['pdf','svg','png']:
        buffer=BytesIO()
        fig.savefig(buffer,format=extension,dpi=600,facecolor='white')
        payload=buffer.getvalue()
        if len(payload)<1000:raise ValueError(f'Incomplete {extension} export for {stem}.')
        destination=figures/(stem+'.'+extension)
        temporary=destination.with_suffix(destination.suffix+'.tmp')
        temporary.write_bytes(payload)
        temporary.replace(destination)
    combined.savefig(fig,facecolor='white')
    plt.close(fig)


def place_circle_labels(ax,x,y,slope,intercept):
    # Same data-coordinate label layout as the current interactive editor.
    dx=np.array([0,-.038044,.042,.040338,-.040338,.042,-.036418,-.021389,.034408,.005852,
                 -.040338,.036418,.034408,0,-.0167,-.021389,-.042,0,-.042,-.031955])
    dy=np.array([-.06,-.025420,0,.016709,.016709,0,.029888,-.051637,.034408,-.059415,
                 .016709,-.029888,.034408,.06,.055053,-.051637,0,-.06,0,-.038937])
    normal=np.array([-slope*np.ptp(x)/np.ptp(y),1.]);normal/=np.linalg.norm(normal)
    offsets=np.column_stack([dx,dy])
    manual={'A':(.0905,.035),'T':(.0960,-.010),'S':(.1020,-.084),
            'Y':(.0968,-.065),'H':(.1340,-.470)}
    positions={}
    for i,(aa,xi,yi) in enumerate(zip(AA,x,y)):
        side=np.sign(yi-(intercept+slope*xi));component=offsets[i]@normal
        if side*component<0:offsets[i]-=2*component*normal
        lab=manual.get(aa,(xi+offsets[i,0]*np.ptp(x),yi+offsets[i,1]*np.ptp(y)))
        ax.annotate(aa,(xi,yi),xytext=lab,ha='center',va='center',fontsize=7,
                    color='white',weight='bold',
                    bbox={'boxstyle':'circle,pad=.2','fc':BLUE,'ec':'#252525','lw':.6},
                    arrowprops={'arrowstyle':'-','color':'#777777','lw':.5},zorder=4)
        positions[aa]=lab
    return positions

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    args=parser.parse_args()
    root=args.root.resolve(); out=root/'outputs'; source=root/'source_data'
    figures=out/'figures'
    for p in [out,source,figures]:p.mkdir(parents=True,exist_ok=True)
    data_path=root/'inputs/residue_RSA_6908_corrected.csv'
    digest=hashlib.sha256(data_path.read_bytes()).hexdigest()
    expected='412039f02ff0fb37730e41eacd1bda6d9defb0922318d927e46088c8a0b147f5'
    if digest!=expected:raise ValueError('The frozen residue dataset checksum does not match.')
    d=pd.read_csv(data_path,keep_default_na=False)
    if set(d.aa)-set(AA):raise ValueError('Unexpected amino-acid codes in frozen input.')
    if d.duplicated(['pdb','chain','residue_number','insertion_code']).any():
        raise ValueError('Duplicate structural residue keys.')
    chain=pd.read_csv(out/'chain_amino_acid_mean_RSA.csv',keep_default_na=False,float_precision='round_trip')
    a=chain.pivot(index=['pdb','chain'],columns='aa',values='mean_RSA').reindex(columns=AA)
    if (len(d),len(a),len(chain))!=(1500829,6342,120822):raise ValueError('Unexpected dataset dimensions.')
    components=pd.read_csv(source/'Table1_full_precision.csv',float_precision='round_trip').set_index('t').loc[AA]
    u=components[['PC1','PC2','PC3']].to_numpy()
    variance=json.loads((source/'pca_variance.json').read_text())
    shares=variance['all_types_percent']
    bins=np.linspace(0,1.5,61); hist=[]
    for aa in AA:
        values=d.loc[d.aa==aa,'rsa'].to_numpy()
        h,_=np.histogram(values,bins=bins)
        for i,n in enumerate(h):hist.append({'aa':aa,'name':NAMES[aa],
            'bin_left':bins[i],'bin_right':bins[i+1],'bin_mid':(bins[i]+bins[i+1])/2,
            'count':int(n),'n_residues_type':len(values),
            'density':float(n/len(values)/(bins[i+1]-bins[i]))})
    h=pd.DataFrame(hist); h.to_csv(source/'Figure1_histogram_bins.csv',index=False)
    statistics={'residue_input_sha256':digest,'n_residues':len(d),'n_chains':len(a),
        'n_chain_type_pairs':len(chain),'RSA_PCA_variance_percent':list(shares[:4]),
        'chain_mean_aggregation':'Original scripts/pool.js sequential sums divided by counts',
        'coordinate_source':'source_data/Table1_full_precision.csv'}
    with PdfPages(figures/'Figures_1-4.pdf') as combined:
        # Figure 1: residue-weighted distributions, not distributions of chain means.
        fig,axs=plt.subplots(4,5,figsize=(7.1,5.5),sharex=True)
        fig.subplots_adjust(left=.085,right=.99,bottom=.09,top=.96,wspace=.40,hspace=.68)
        for ax,aa in zip(axs.flat,AA):
            hh=h[h.aa==aa]
            ax.stairs(hh.density.to_numpy(),bins,fill=True,color=BLUE,lw=.45)
            ax.set_title(f'{aa}  {NAMES[aa]}',pad=3,fontweight='bold',fontsize=8)
            ax.set_xlim(0,1.5);ax.set_ylim(bottom=0);ax.set_xticks([0,.5,1,1.5])
            ax.tick_params(length=2,pad=2)
        fig.supxlabel('Residue relative solvent accessibility (RSA)',y=.015,fontsize=9)
        fig.supylabel('Probability density',x=.012,fontsize=9)
        save(fig,'Figure1_residue_RSA',figures,combined)
        # Figure 2: full-precision PCA coordinates and sample SD across chains.
        sd=a.std(ddof=1).to_numpy(); pearson=stats.pearsonr(sd,u[:,1]);spearman=stats.spearmanr(sd,u[:,1])
        fit=stats.linregress(sd,u[:,1])
        f2=pd.DataFrame({'aa':AA,'n_chains_with_aa':a.count().to_numpy(),
            'sd_chain_mean_RSA':sd,'PC2':u[:,1]})
        fig,ax=plt.subplots(figsize=(4.2,3.7));fig.subplots_adjust(left=.17,bottom=.16,right=.97,top=.96)
        ax.scatter(sd,u[:,1],s=17,color='#252525',zorder=3)
        xx=np.linspace(sd.min()-.001,sd.max()+.001,200)
        ax.plot(xx,fit.intercept+fit.slope*xx,color=ORANGE,lw=1.5,zorder=1)
        ax.set(xlabel='SD of per-chain mean RSA',ylabel='PC2 coordinate',
               xlim=(.062,.158),ylim=(-.67,.43))
        positions=place_circle_labels(ax,sd,u[:,1],fit.slope,fit.intercept)
        f2['label_x']=[positions[t][0] for t in AA];f2['label_y']=[positions[t][1] for t in AA]
        ax.text(.98,.98,f'Pearson r = {pearson.statistic:.3f}',
            transform=ax.transAxes,ha='right',va='top',fontsize=8,fontweight='bold')
        f2.to_csv(source/'Figure2_PC2_dispersion.csv',index=False)
        pd.DataFrame({'sd_chain_mean_RSA':xx,'fitted_PC2':fit.intercept+fit.slope*xx}).to_csv(source/'Figure2_fit_line.csv',index=False)
        statistics['Figure2']={'Pearson_r':float(pearson.statistic),'Pearson_P_two_sided':float(pearson.pvalue),
            'Spearman_rho':float(spearman.statistic),'Spearman_P_two_sided':float(spearman.pvalue),
            'n_amino_acid_types':20,'line_slope':float(fit.slope),'line_intercept':float(fit.intercept)}
        save(fig,'Figure2_PC2_dispersion',figures,combined)
        # Figure 3: exactly zero plus ten positive intervals; the final bin is >90%.
        f3=pd.read_csv(source/'TableS1_full_precision.csv',float_precision='round_trip')
        f3=f3.rename(columns={'interval':'interval_label','pairs':'n_chain_type_pairs','pct':'percentage_all_pairs','cv':'CV_across_types'})
        f3.insert(0,'interval_order',range(11))
        f3.to_csv(source/'Figure3_divergence.csv',index=False)
        pd.read_csv(source/'TableS1_fractions_by_type.csv',float_precision='round_trip').to_csv(source/'Figure3_interval_fractions_by_type.csv',index=False)
        fig,ax=plt.subplots(figsize=(7.1,3.3));fig.subplots_adjust(left=.10,bottom=.30,right=.99,top=.94)
        # Source CSV retains raw CV; the displayed axis is 100 x CV.
        ax.plot(range(11),100*f3.CV_across_types,color=BLUE,lw=1.5)
        ax.scatter(range(11),100*f3.CV_across_types,s=30,color=ORANGE,
                   edgecolor='#252525',linewidth=.5,zorder=3)
        labels=['0']+[f'{i/10:g}–{(i+1)/10:g}' for i in range(9)]+['>0.9']
        ax.set(xticks=range(11),xticklabels=labels,xlabel='Per-chain mean RSA interval',
            ylabel='Coefficient of variation (%)',ylim=(0,1.1*100*f3.CV_across_types.max()))
        ax.tick_params(axis='x',labelrotation=45)
        for tick in ax.get_xticklabels():tick.set_ha('right')
        save(fig,'Figure3_divergence',figures,combined)
        # Figure 4: both kingdoms are in one joint component space.
        assignment=json.loads((out/'assign.json').read_text())
        labels=np.array([assignment[(p+c).lower()]['k'] for p,c in a.index])
        euk=a.loc[labels=='euk'];prok=a.loc[labels=='prok']
        if (len(euk),len(prok))!=(1778,2877):raise ValueError('Unexpected kingdom counts.')
        joint_table=pd.read_csv(source/'TableS2_full_precision.csv',float_precision='round_trip').set_index('t').loc[AA]
        joint=np.vstack([joint_table[['PC1_euk','PC2_euk']].to_numpy(),joint_table[['PC1_prok','PC2_prok']].to_numpy()])
        joint_shares=variance['joint_kingdom_percent']
        f4=pd.DataFrame({'aa':AA})
        fig,axs=plt.subplots(1,2,figsize=(7.1,5.2),sharey=True)
        fig.subplots_adjust(left=.09,right=.98,bottom=.14,top=.89,wspace=.36)
        order=np.arange(19,-1,-1)
        for k,ax in enumerate(axs):
            x=joint[20:,k];y=joint[:20,k];yy=np.arange(20)
            ax.hlines(yy,x[order],y[order],color='#777777',lw=.8)
            ax.scatter(y[order],yy,s=25,color=BLUE,edgecolor='#252525',lw=.6,
                       label='Eukaryotes',zorder=3)
            ax.scatter(x[order],yy,s=25,facecolor='white',edgecolor=ORANGE,lw=1,
                       label='Prokaryotes',zorder=3)
            ax.axvline(0,color='#aaaaaa',ls=':',lw=.6,zorder=0)
            ax.set(yticks=yy,yticklabels=np.array(AA)[order],ylim=(-.7,19.7),
                   xlabel=f'PC{k+1} coordinate\nVariance explained: {joint_shares[k]:.2f}%')
            ax.set_title('ab'[k],loc='left',fontweight='bold')
            f4[f'PC{k+1}_prok']=x;f4[f'PC{k+1}_euk']=y
            f4[f'delta_PC{k+1}_euk_minus_prok']=y-x
        axs[0].set_ylabel('Amino-acid type')
        axs[1].legend(loc='lower right',bbox_to_anchor=(1,1.02),ncol=2,
                      frameon=False,fontsize=7)
        f4.to_csv(source/'Figure4_kingdom_coordinates.csv',index=False)
        f4.to_csv(source/'TableS2_joint_kingdom_coordinates.csv',index=False)
        statistics['Figure4']={'n_eukaryotic_chains':len(euk),'n_prokaryotic_chains':len(prok),
            'joint_variance_percent':list(joint_shares[:3]),
            'mean_delta_PC1_euk_minus_prok':float((joint[:20,0]-joint[20:,0]).mean()),
            'mean_delta_PC2_euk_minus_prok':float((joint[:20,1]-joint[20:,1]).mean())}
        save(fig,'Figure4_kingdom_coordinates',figures,combined)
    statistics['runtime']={'Python':sys.version.split()[0],'NumPy':np.__version__,'pandas':pd.__version__,
        'SciPy':scipy.__version__,'Matplotlib':matplotlib.__version__,'platform':platform.platform()}
    (source/'figure_statistics.json').write_text(json.dumps(statistics,indent=2)+'\n')
    print('Created four PDF/SVG/600-dpi PNG figures, a combined PDF and source CSVs.')

if __name__=='__main__':main()
