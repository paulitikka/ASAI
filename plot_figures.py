"""Regenerate Figures 1-4 and their source data from the frozen RSA input.

Run: python plot_figures.py
Only local, bundled files are read. No manuscript values are fitted or forced.
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
BLUE = '#477B8A'
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
    g=d.groupby(['pdb','chain','aa'],sort=True).rsa.agg(['mean','count'])
    a=g['mean'].unstack('aa').reindex(columns=AA)
    counts=g['count'].unstack('aa').reindex(columns=AA)
    if (len(d),len(a),len(g))!=(1500829,6342,120822):raise ValueError('Unexpected dataset dimensions.')
    # Use the original sequential pooling for every chain-mean analysis. This
    # preserves the numerical convention at exact interval boundaries; compensated
    # pandas summation may otherwise move a handful of pairs across a boundary.
    pool=json.loads((root/'work/pooled-acc.json').read_text())
    pool_order='ARNDCQEGHILKMFPSTWYV'
    canonical=[]
    for pdb,chain in a.index:
        item=pool[(pdb+chain).lower()]
        canonical.append([item['sum'][pool_order.index(t)]/item['cnt'][pool_order.index(t)]
            if item['cnt'][pool_order.index(t)] else np.nan for t in AA])
    canonical=np.asarray(canonical)
    if np.nanmax(np.abs(canonical-a.to_numpy()))>1e-12:
        raise ValueError('Original pooling and independent aggregation disagree.')
    a=pd.DataFrame(canonical,index=a.index,columns=AA)
    g['mean']=a.stack().reindex(g.index)
    g.rename(columns={'mean':'mean_RSA','count':'n_residues'}).to_csv(out/'chain_amino_acid_mean_RSA.csv')
    q=np.nanquantile(a.to_numpy(),np.arange(1,100)/100,axis=0).T
    u,shares=decompose(q)
    mean=a.mean().to_numpy()
    if np.dot(u[:,0],mean)<0:u[:,0]*=-1
    if u[AA.index('K'),1]<u[AA.index('G'),1]:u[:,1]*=-1
    archived=json.loads((out/'decompose.json').read_text())
    ref=pd.DataFrame(archived['table1']).set_index('t').loc[AA]
    if np.dot(u[:,2],ref.PC3.to_numpy())<0:u[:,2]*=-1
    table1=pd.DataFrame({'aa':AA,'n_chains_with_aa':a.count().to_numpy(),
        'n_residues':d.groupby('aa').size().reindex(AA).to_numpy(),
        'mean_RSA_equal_chain_weights':mean,
        'mean_RSA_equal_residue_weights':d.groupby('aa').rsa.mean().reindex(AA).to_numpy(),
        'PC1':u[:,0],'PC2':u[:,1],'PC3':u[:,2]})
    table1.to_csv(source/'Table1_RSA_components.csv',index=False)
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
        'n_chain_type_pairs':len(g),'RSA_PCA_variance_percent':shares[:4].tolist(),
        'chain_mean_aggregation':'Original scripts/pool.js sequential sums divided by counts',
        'PC2_max_absolute_difference_from_archived_rounded_coordinates':
            float(np.max(np.abs(u[:,1]-ref.PC2.to_numpy())))}
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
        ax.scatter(sd,u[:,1],s=25,color=BLUE,edgecolor='white',linewidth=.45,zorder=3)
        xx=np.linspace(sd.min(),sd.max(),200)
        ax.plot(xx,fit.intercept+fit.slope*xx,color='#777777',lw=.8,zorder=1)
        ax.set(xlabel='SD of chain-specific mean RSA',ylabel='PC2 coordinate')
        ax.margins(x=.17,y=.16)
        positions=place_labels(ax,sd,u[:,1],AA)
        f2['label_x']=[positions[t][0] for t in AA];f2['label_y']=[positions[t][1] for t in AA]
        exponent=int(np.floor(np.log10(pearson.pvalue)));mantissa=pearson.pvalue/10**exponent
        ax.text(.98,.98,f'Pearson r = {pearson.statistic:.3f}\n'
            +rf'$P = {mantissa:.2f}\times 10^{{{exponent}}}$'+'\n'+f'n = {len(AA)} amino-acid types',
            transform=ax.transAxes,ha='right',va='top',fontsize=7.4)
        f2.to_csv(source/'Figure2_PC2_dispersion.csv',index=False)
        pd.DataFrame({'sd_chain_mean_RSA':xx,'fitted_PC2':fit.intercept+fit.slope*xx}).to_csv(source/'Figure2_fit_line.csv',index=False)
        statistics['Figure2']={'Pearson_r':float(pearson.statistic),'Pearson_P_two_sided':float(pearson.pvalue),
            'Spearman_rho':float(spearman.statistic),'Spearman_P_two_sided':float(spearman.pvalue),
            'n_amino_acid_types':20,'line_slope':float(fit.slope),'line_intercept':float(fit.intercept)}
        save(fig,'Figure2_PC2_dispersion',figures,combined)
        # Figure 3: exactly zero plus ten positive intervals; the final bin is >90%.
        divergence=[]; fractions=[]; v=a.to_numpy();valid=np.isfinite(v).sum(axis=0)
        definitions=[(0,0,'0')]+[(i/10,(i+1)/10 if i<9 else np.inf,
            f'({i*10}, {(i+1)*10}]' if i<9 else '>90') for i in range(10)]
        for i,(lo,hi,label) in enumerate(definitions):
            mask=(v==0) if i==0 else ((v>lo)&(v<=hi))
            bin_counts=mask.sum(axis=0);fraction=bin_counts/valid
            cv=float(fraction.std(ddof=1)/fraction.mean());n=int(bin_counts.sum())
            divergence.append({'interval_order':i,'interval_label':label,'n_chain_type_pairs':n,
                'percentage_all_pairs':100*n/len(g),'CV_across_types':cv})
            for j,aa in enumerate(AA):fractions.append({'interval_order':i,'interval_label':label,
                'aa':aa,'pairs_in_interval':int(bin_counts[j]),'chains_with_type':int(valid[j]),
                'fraction_within_type':float(fraction[j])})
        f3=pd.DataFrame(divergence)
        f3.to_csv(source/'Figure3_divergence.csv',index=False)
        pd.DataFrame(fractions).to_csv(source/'Figure3_interval_fractions_by_type.csv',index=False)
        fig,ax=plt.subplots(figsize=(7.1,3.3));fig.subplots_adjust(left=.10,bottom=.30,right=.99,top=.94)
        ax.plot(range(11),f3.CV_across_types,'o-',color=BLUE,lw=1.1,ms=4)
        ax.set(xticks=range(11),xticklabels=f3.interval_label,xlabel='Chain-specific mean RSA interval (%)',
            ylabel='Coefficient of variation across types',ylim=(0,3.12))
        ax.tick_params(axis='x',labelrotation=45)
        for tick in ax.get_xticklabels():tick.set_ha('right')
        save(fig,'Figure3_divergence',figures,combined)
        # Figure 4: both kingdoms are in one joint component space.
        assignment=json.loads((out/'assign.json').read_text())
        labels=np.array([assignment[(p+c).lower()]['k'] for p,c in a.index])
        euk=a.loc[labels=='euk'];prok=a.loc[labels=='prok']
        if (len(euk),len(prok))!=(1778,2877):raise ValueError('Unexpected kingdom counts.')
        rows=np.vstack([np.nanquantile(euk,np.arange(1,100)/100,axis=0).T,
                        np.nanquantile(prok,np.arange(1,100)/100,axis=0).T])
        joint,joint_shares=decompose(rows)
        if np.dot(joint[:,0],np.r_[euk.mean().to_numpy(),prok.mean().to_numpy()])<0:joint[:,0]*=-1
        if (joint[AA.index('K'),1]+joint[20+AA.index('K'),1])<(joint[AA.index('G'),1]+joint[20+AA.index('G'),1]):joint[:,1]*=-1
        f4=pd.DataFrame({'aa':AA})
        fig,axs=plt.subplots(1,2,figsize=(7.1,3.6));fig.subplots_adjust(left=.09,right=.98,bottom=.20,top=.89,wspace=.36)
        for k,ax in enumerate(axs):
            x=joint[20:,k];y=joint[:20,k];lo=min(x.min(),y.min());hi=max(x.max(),y.max());pad=(hi-lo)*.18
            limits=(lo-pad,hi+pad);ax.set(xlim=limits,ylim=limits,aspect='equal',
                xlabel=f'Prokaryotic PC{k+1} coordinate',ylabel=f'Eukaryotic PC{k+1} coordinate')
            ax.plot(limits,limits,color='#999999',ls='--',lw=.7,zorder=1)
            ax.scatter(x,y,s=24,color=BLUE,edgecolor='white',linewidth=.45,zorder=3)
            ax.set_title(f'{chr(97+k)}  PC{k+1}',loc='left',fontweight='bold',pad=8)
            positions=place_labels(ax,x,y,AA)
            f4[f'PC{k+1}_prok']=x;f4[f'PC{k+1}_euk']=y;f4[f'delta_PC{k+1}_euk_minus_prok']=y-x
            f4[f'label_x_PC{k+1}']=[positions[t][0] for t in AA]
            f4[f'label_y_PC{k+1}']=[positions[t][1] for t in AA]
        f4.to_csv(source/'Figure4_kingdom_coordinates.csv',index=False)
        f4.to_csv(source/'TableS2_joint_kingdom_coordinates.csv',index=False)
        statistics['Figure4']={'n_eukaryotic_chains':len(euk),'n_prokaryotic_chains':len(prok),
            'joint_variance_percent':joint_shares[:3].tolist(),
            'mean_delta_PC1_euk_minus_prok':float((joint[:20,0]-joint[20:,0]).mean()),
            'mean_delta_PC2_euk_minus_prok':float((joint[:20,1]-joint[20:,1]).mean())}
        save(fig,'Figure4_kingdom_coordinates',figures,combined)
    # A readable annotation table; original input, including its legacy taxid field, is unchanged.
    annotation=pd.read_csv(root/'inputs/selected_6908_organisms.csv',keep_default_na=False)
    annotation=annotation.rename(columns={'taxid':'legacy_taxid_annotation'})
    annotation['classification_key']=(annotation.pdb+annotation.chain).str.lower()
    annotation['kingdom_classification']=annotation.classification_key.map(lambda k:assignment.get(k,{}).get('k','unknown'))
    structural_counts=d.groupby(['pdb','chain']).size()
    lookup={(str(p).lower(),str(c)):int(n) for (p,c),n in structural_counts.items()}
    annotation['n_RSA_residues']=[lookup.get((str(p).lower(),str(c)),0) for p,c in zip(annotation.pdb,annotation.chain)]
    annotation['has_processed_RSA']=annotation.n_RSA_residues>0
    annotation.to_csv(out/'chain_selection_organism_classification.csv',index=False)
    models=json.loads((out/'admissible_parameters.json').read_text())
    pd.DataFrame(models['parameters']).to_csv(source/'Table2_model_parameters.csv',index=False)
    model_labels={'length':'Length alone','composition':'Length + composition',
        'charge':'Length + net charge + signed charge balance',
        'hydrophobicity':'Length + sequence-encoded hydrophobicity','all':'All five parameters'}
    pd.DataFrame([{'model_key':key,'model':model_labels[key],
        'mean_held_out_R2':value['meanHeldOutR2'],'unexplained_percent':value['unexplainedPct'],
        'n_type_split_evaluations':value['n']} for key,value in models['models'].items()]
        ).to_csv(source/'Table3_joint_models.csv',index=False)
    pd.DataFrame(json.loads((out/'s6_windows.json').read_text())).to_csv(
        source/'TableS3_length_restricted_comparison.csv',index=False)
    statistics['runtime']={'Python':sys.version.split()[0],'NumPy':np.__version__,'pandas':pd.__version__,
        'SciPy':scipy.__version__,'Matplotlib':matplotlib.__version__,'platform':platform.platform()}
    (source/'figure_statistics.json').write_text(json.dumps(statistics,indent=2)+'\n')
    print(json.dumps(statistics,indent=2))
    print('Created four PDF/SVG/600-dpi PNG figures, a combined PDF and source CSVs.')

if __name__=='__main__':main()
