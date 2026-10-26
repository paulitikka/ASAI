from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.stats import beta,norm
P=Path(__file__).resolve().parent;root=P.parent;d=np.load(P/'rsa_distribution_data.npz');a=d['a'];n=d['n'];u=d['u'];v=d['v'];q=d['q'];s=d['s'];rng=np.random.default_rng(20261005)
def cc(x,y):return float(np.corrcoef(x,y)[0,1])
def quant(x):return np.nanquantile(x,np.arange(1,100)/100,axis=0).T
def pc(x):return np.linalg.svd(x-x.mean(0),full_matrices=False)
def residual(y,U):
 X=np.column_stack([np.ones(len(y)),U]);return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
R={};R['variants']={}
variants={f'min_residues_per_chain_type_{k}':np.where(n>=k,a,np.nan) for k in [3,5,10,20]};variants['all_20_types_present']=a[np.isfinite(a).all(1)]
length=np.nansum(n,axis=1)
for k in [100,200,300]:variants[f'chain_length_at_least_{k}']=a[length>=k]
for name,x in variants.items():
 qq=quant(x);uu,ss,vv=pc(qq);rr=cc(uu[:,2],u[:,2]);sgn=1 if rr>=0 else -1
 R['variants'][name]={'min_chains_per_AA':int(np.isfinite(x).sum(0).min()),'PC3_correlation':abs(rr),'PC3_share':float(100*ss[2]**2/(ss@ss)),'loading_correlation':cc(vv[2]*sgn,v[2])}
for low,high in [(5,95),(10,90)]:
 uu,ss,vv=pc(q[:,low-1:high]);R['variants'][f'quantiles_{low}_{high}']={'PC3_correlation':abs(cc(uu[:,2],u[:,2])),'PC3_share':float(100*ss[2]**2/(ss@ss))}
# Moment-matched two-parameter beta benchmark, same support convention as earlier work.
mean=np.nanmean(a,0);sd=np.nanstd(a,0,ddof=1);m=mean/1.5;var=(sd/1.5)**2;t=m*(1-m)/var-1
bq=beta.ppf(np.arange(1,100)[None,:]/100,(m*t)[:,None],((1-m)*t)[:,None])*1.5
bu,bs,bv=pc(bq);R['two_parameter_beta']={'PC3_correlation':abs(cc(bu[:,2],u[:,2])),'PC3_share':float(100*bs[2]**2/(bs@bs)),'observed_PC3_vs_beta_asymmetry_partial':cc(u[:,2],residual(bq[:,94]+bq[:,4]-2*bq[:,49],u[:,:2])),'note':'Descriptive counterexample; not a validated generative model; support [0,1.5].'}
# Independent chain halves: association of training PC3 with held-out asymmetry after removing training PC1/2.
cv=[]
for b in range(100):
 ix=rng.permutation(len(a));qa=quant(a[ix[:len(ix)//2]]);qb=quant(a[ix[len(ix)//2:]]);uu,ss,vv=pc(qa)
 if cc(uu[:,2],u[:,2])<0:uu[:,2]*=-1
 asym=qb[:,94]+qb[:,4]-2*qb[:,49];cv.append(cc(uu[:,2],residual(asym,uu[:,:2])))
R['split_chain_validation']={'replicates':100,'partial_asymmetry_correlation_quantiles_025_50_975':np.quantile(cv,[.025,.5,.975]).tolist()}
# Exploratory property scan, control the maximum raw correlation within each searched family.
records=json.load(open(root/'inputs/aaindex_records.json'));audit=pd.read_csv(root/'intrinsic/AAindex_553_audit.csv').set_index('accession');mat=np.column_stack([x['values'] for x in records]);mat=mat-mat.mean(0);mat=mat/np.linalg.norm(mat,axis=0);y=u[:,2];obs=y@mat
perms=np.array([rng.permutation(20) for _ in range(20000)]);null=y[perms]@mat
rows=[];sel=np.array([bool(audit.loc[x['accession'],'selected_primary']) for x in records]);maxall=np.max(abs(null),axis=1);maxsel=np.max(abs(null[:,sel]),axis=1)
for i,x in enumerate(records):
 c=x['accession'];raw=np.asarray(x['values']);rows.append({'accession':c,'description':audit.loc[c,'description'],'selected_intrinsic':bool(sel[i]),'r_PC3':float(obs[i]),'partial_r_PC3':cc(y,residual(raw,u[:,:2])),'p_max_553':float((1+(maxall>=abs(obs[i])).sum())/20001),'p_max_51':float((1+(maxsel>=abs(obs[i])).sum())/20001) if sel[i] else None})
pd.DataFrame(rows).sort_values('r_PC3',key=abs,ascending=False).to_csv(P/'AAindex_PC3_correlations.csv',index=False)
R['best_intrinsic']=max([r for r in rows if r['selected_intrinsic']],key=lambda r:abs(r['r_PC3']));R['best_all']=max(rows,key=lambda r:abs(r['r_PC3']))
R['normal_quantile_curvature_correlation']=cc(v[2],norm.ppf(np.arange(1,100)/100)**2-1)
(P/'validation_results.json').write_text(json.dumps(R,indent=2));print(json.dumps(R,indent=2))
