"""Independent raw-data, PCA and linear-model checks for the r3 release.

The interval benchmark reconstructs the integer/half-integer ASA lattice and
uses integer comparisons, independently of the production floating tolerance.
The models are independently solved with NumPy using the archived split ranks.
"""
from pathlib import Path
import json,math
import numpy as np
import pandas as pd

P=Path(__file__).resolve().parent
AA=list('ACDEFGHIKLMNPQRSTVWY');BASE=list('ARNDCQEGHILKMFPSTWYV')
MAX={'A':129,'R':274,'N':195,'D':193,'C':167,'Q':225,'E':223,'G':104,'H':224,'I':197,'L':201,'K':236,'M':224,'F':240,'P':159,'S':155,'T':172,'W':285,'Y':263,'V':174}

def main():
    raw=pd.read_csv(P/'inputs/residue_RSA_6908_corrected.csv',keep_default_na=False,float_precision='round_trip')
    assert len(raw)==1500829 and not raw.duplicated(['pdb','chain','residue_number','insertion_code']).any()
    sizes=raw.groupby(['pdb','chain']).size();assert len(sizes)==6342
    maximum=raw.aa.map(MAX).to_numpy();raw['twice_ASA']=np.rint(2*maximum*raw.rsa).astype(np.int64)
    lattice_error=float(np.max(abs(raw.rsa-raw.twice_ASA/(2*maximum))))
    assert lattice_error<3e-16
    g=raw.groupby(['pdb','chain','aa'],sort=True).agg(sum_ASA=('twice_ASA','sum'),n=('rsa','size'),mean_RSA=('rsa','mean'))
    denominators=np.array([2*MAX[t]*int(n) for t,n in zip(g.index.get_level_values('aa'),g.n)],dtype=np.int64)
    # Exact nominal intervals for the documented integer-ASA normalization.
    integer_bins=np.minimum(10,(10*g.sum_ASA.to_numpy()+denominators-1)//denominators)
    counts=pd.crosstab(integer_bins,g.index.get_level_values('aa')).reindex(index=range(11),columns=AA,fill_value=0)
    prod=json.loads((P/'outputs/decompose.json').read_text())
    assert np.array_equal(counts.to_numpy(),np.array(prod['intervalCountsByType']))
    fr=counts/counts.sum(0);cv=fr.std(1,ddof=1)/fr.mean(1)
    assert np.allclose(cv,[r['cv'] for r in prod['S5FullPrecision']],atol=1e-12,rtol=0)
    a=g.mean_RSA.unstack('aa').reindex(columns=AA)
    canonical=pd.read_csv(P/'outputs/chain_amino_acid_mean_RSA.csv',float_precision='round_trip')
    archived=canonical.pivot(index=['pdb','chain'],columns='aa',values='mean_RSA').reindex(columns=AA)
    aggregation_error=float(np.nanmax(abs(a.to_numpy()-archived.to_numpy())))
    assert aggregation_error<1e-12
    q=np.nanquantile(a,np.arange(1,100)/100,axis=0).T
    u,s,v=np.linalg.svd(q-q.mean(0),full_matrices=False)
    ref=pd.DataFrame(prod['table1FullPrecision']).set_index('t').loc[AA]
    for k in range(3):
        if np.dot(u[:,k],ref[f'PC{k+1}'])<0:u[:,k]*=-1
    pca_error=float(np.max(abs(u[:,:3]-ref[['PC1','PC2','PC3']].to_numpy())))
    assert pca_error<5e-9
    assert np.allclose(100*s*s/(s@s),prod['shares1FullPrecision'],atol=1e-8)

    pool=json.loads((P/'work/pooled-acc.json').read_text());predictors=json.loads((P/'outputs/chain_predictors.json').read_text())
    selected=pd.read_csv(P/'inputs/selected_6908_organisms.csv',keep_default_na=False)
    sequences={(r.pdb+r.chain).lower():r.sequence for r in selected.itertuples()}
    KD={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-.8,'T':-.7,'W':-.9,'Y':-1.3,'V':4.2}
    independent={}
    for key,r in predictors.items():
        c=pool[key];n=np.array(c['cnt'],float);L=c['L'];seq=sequences[key]
        h=np.array([KD.get(t,0) for t in seq]);blocks=[seq[i:i+20] for i in range(0,len(seq)-19,20)]
        charged={t:n[BASE.index(t)] for t in ['D','E','K','R','H']}
        fresh={'length':math.log10(L),'comp':((n+.5)/(L+10)).tolist(),
               'charge':(charged['K']+charged['R']+.1*charged['H']-charged['D']-charged['E'])/L,
               'asymmetry':(charged['K']+charged['R']+charged['H']-charged['D']-charged['E'])/L,
               'hydroMean':float(h.mean()),'hydroWindow':float(np.convolve(h,np.ones(11)/11,'valid').max()) if len(h)>=11 else 0,
               'lowComplexity':sum(len(set(s))<=12 for s in blocks)/len(blocks) if blocks else 0}
        for field in fresh:assert np.allclose(fresh[field],r[field],atol=1e-11,rtol=0),(key,field)
        independent[key]=fresh
    design=pd.read_csv(P/'design/train_test_splits.csv',dtype={'chain_key':str,'aa':str})
    diagnostics=pd.read_csv(P/'source_data/Table3_fit_diagnostics.csv').set_index(['split','aa','model'])
    max_error=0.;computed=[]
    def features(r,name):
        length=[r['length']];comp=r['comp'];charge=[r['charge'],r['asymmetry']];hydro=[r['hydroMean'],r['hydroWindow'],r['lowComplexity']]
        body={'length':length,'composition':comp+length,'charge':charge+length,'hydrophobicity':hydro+length,'all':comp+charge+hydro+length}[name]
        return [1.]+body
    for aa in BASE:
        subset=design[design.aa==aa].set_index('chain_key')
        for split in range(1,21):
            column=f'permutation_rank_{split:02d}';order=subset.sort_values(column).index.tolist()
            assert sorted(subset[column])==list(range(1,len(order)+1))
            half=len(order)//2;train,test=order[:half],order[half:]
            assert not set(train)&set(test)
            response=lambda keys:np.array([pool[k]['sum'][BASE.index(aa)]/pool[k]['cnt'][BASE.index(aa)] for k in keys])
            ytr,yte=response(train),response(test)
            for name in ['length','composition','charge','hydrophobicity','all']:
                X=np.array([features(independent[k],name) for k in train]);Z=np.array([features(independent[k],name) for k in test])
                A=X.T@X;lam=1e-9*np.trace(A)/A.shape[0]
                beta=np.linalg.solve(A+lam*np.eye(len(A)),X.T@ytr)
                r2=1-np.sum((yte-Z@beta)**2)/np.sum((yte-yte.mean())**2)
                original=float(diagnostics.loc[(split,aa,name),'R2']);max_error=max(max_error,abs(float(r2)-original))
                computed.append({'split':split,'aa':aa,'model':name,'independent_R2':r2,'JavaScript_R2':original})
    assert max_error<1e-6,max_error
    calculations=pd.DataFrame(computed)
    for name,x in calculations.groupby('model'):
        assert f'{x.independent_R2.mean():.4f}'==f'{x.JavaScript_R2.mean():.4f}',name
    calculations.to_csv(P/'outputs/independent_model_fit_comparison.csv',index=False)
    report={'status':'PASS','version':'2026.10.26-r3','checks':['raw residue keys and dimensions',
        'integer-ASA interval assignment independent of floating tolerance','independent interval fractions and CV',
        'independent compensated aggregation and SVD','all 6287 predictor records independently recalculated',
        'all 2000 held-out fits independently solved with archived split ranks'],
        'RSA_lattice_max_error':lattice_error,'aggregation_max_error':aggregation_error,
        'PCA_unit_coordinate_max_error':pca_error,'independent_fit_R2_max_error':max_error,
        'n_independent_model_fits':len(computed)}
    (P/'outputs/independent_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
