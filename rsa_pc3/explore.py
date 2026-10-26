from pathlib import Path
import pandas as pd,numpy as np
import hashlib
from scipy.stats import skew,kurtosis
P=Path(__file__).resolve().parent;root=P.parent
assert hashlib.sha256((P/'input/residue_RSA_6908_corrected.csv').read_bytes()).hexdigest()=='412039f02ff0fb37730e41eacd1bda6d9defb0922318d927e46088c8a0b147f5'
order=list('ARNDCQEGHILKMFPSTWYV');d=pd.read_csv(P/'input/residue_RSA_6908_corrected.csv');g=d.groupby(['pdb','chain','aa']).rsa.agg(['mean','count']);a=g['mean'].unstack('aa').reindex(columns=order);n=g['count'].unstack('aa').reindex(columns=order)
q=np.nanquantile(a,np.arange(1,100)/100,axis=0).T;xc=q-q.mean(0);u,s,v=np.linalg.svd(xc,full_matrices=False)
ref=pd.read_csv(P/'rsa_coordinates.csv',index_col=0).loc[order]
for k in range(3):
 if np.corrcoef(u[:,k],ref[f'RSA_PC{k+1}'])[0,1]<0:u[:,k]*=-1;v[k]*=-1
m=pd.DataFrame(index=order)
m['mean']=a.mean();m['sd']=a.std();m['skewness']=skew(a,axis=0,nan_policy='omit',bias=False);m['excess_kurtosis']=kurtosis(a,axis=0,nan_policy='omit',bias=False)
m['median']=q[:,49];m['mean_minus_median']=m['mean']-m['median'];m['Bowley_skew']=(q[:,74]+q[:,24]-2*q[:,49])/(q[:,74]-q[:,24]);m['tail_asymmetry90']=(q[:,94]+q[:,4]-2*q[:,49])/(q[:,94]-q[:,4]);m['q95_plus_q05_minus_2median']=q[:,94]+q[:,4]-2*q[:,49];m['q95_minus_median']=q[:,94]-q[:,49];m['median_minus_q05']=q[:,49]-q[:,4];m['tail_weight']=(q[:,94]-q[:,4])/(q[:,74]-q[:,24]);m['chains_with_aa']=a.count();m['mean_count_per_chain']=n.mean();m['fraction_singletons']=(n==1).sum()/n.count();m['fraction_zero_chainmean']=(a==0).sum()/a.count()
controls=np.column_stack([np.ones(20),u[:,:2]]);out=[]
for c in m:
 y=m[c].values;r=y-controls@np.linalg.lstsq(controls,y,rcond=None)[0];out.append([c,np.corrcoef(y,u[:,2])[0,1],np.corrcoef(r,u[:,2])[0,1]])
print(pd.DataFrame(out,columns=['property','r_PC3','partial_PC3']).round(4).to_string(index=False));print('LOADINGS',[(i+1,round(v[2,i],4)) for i in [0,4,9,24,49,74,89,94,98]])
print(m[['skewness','Bowley_skew','tail_asymmetry90','mean_count_per_chain']].round(3).to_string());m['RSA_PC3']=u[:,2];m.to_csv(P/'distribution_descriptors.csv');np.savez(P/'rsa_distribution_data.npz',a=a.values,n=n.values,q=q,u=u,s=s,v=v,order=order)
