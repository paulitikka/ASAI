import json
from pathlib import Path
import numpy as np,pandas as pd
P=Path(__file__).resolve().parent;rng=np.random.default_rng(10604)
X=pd.read_csv(P/'selected_intrinsic_values.csv',index_col=0);Z=(X-X.mean())/X.std(ddof=1)
def pa(x,B=1500):
 x=np.asarray(x);real=np.linalg.eigvalsh(x@x.T)[::-1];s=[]
 for b in range(B):
  z=np.take_along_axis(x,np.argsort(rng.random(x.shape),axis=0),axis=0);s.append(np.linalg.eigvalsh(z@z.T)[::-1])
 q=np.quantile(s,.95,axis=0);return {'retained':next((i for i in range(19) if real[i]<=q[i]),19),'PC3_p':float((1+(np.array(s)[:,2]>=real[2]).sum())/(B+1)),'spectrum_percent':(100*real/real.sum()).tolist(),'PC3_percent':float(100*real[2]/real.sum()),'PC3_null95_percent':float(100*q[2]/real.sum())}
fam=np.array([c[:4] for c in Z.columns]);out={'description':'Post-primary diagnostic: leave one accession-author family out; 1500 permutations each. Families proxy source redundancy, not independent studies.','families':{}}
for f in sorted(set(fam)):out['families'][f]=pa(Z.loc[:,fam!=f].values)
W=Z.values/np.sqrt(np.array([sum(fam==f) for f in fam]));out['equal_family_weight']=pa(W,5000)
(P/'family_sensitivity.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
