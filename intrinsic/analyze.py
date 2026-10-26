from audit import *
rng=np.random.default_rng(20261005); B=5000
D={r['accession']:np.array(r['values']) for r in R}
def centered(x):return x-x.mean(0)
def decomp(x):return np.linalg.svd(centered(x),full_matrices=False)
def eig(x):return np.maximum(np.linalg.eigvalsh(x@x.T)[::-1],0)
def pa(x,groups=None):
 real=eig(x);sim=[]
 for b in range(B):
  if groups is None:z=np.take_along_axis(x,np.argsort(rng.random(x.shape),axis=0),axis=0)
  else:
   z=x.copy()
   for g in np.unique(groups):z[:,groups==g]=x[rng.permutation(20)][:,groups==g]
  sim.append(eig(z))
 sim=np.array(sim);q=np.quantile(sim,.95,axis=0);keep=next((i for i in range(19) if real[i]<=q[i]),19)
 return dict(retained=keep,observed_percent=real/real.sum()*100,null95_percent=q/real.sum()*100,p=(1+(sim>=real).sum(0))/(B+1))
A=pd.read_csv(P/'rsa_coordinates.csv').set_index('aa').loc[order,[f'RSA_PC{i}' for i in range(1,4)]].values
perms=np.array([rng.permutation(20) for _ in range(20000)])
def cross(I):
 G=I.T@A;l,c,r=np.linalg.svd(G,full_matrices=True);U=I@l;V=A@r.T;ps=[]
 for k in range(len(c)):
  vals=np.linalg.svd(np.einsum('ni,bnj->bij',U[:,k:],V[perms,k:]),compute_uv=False)
  stats=-np.log1p(-np.minimum(vals**2,1-1e-14)).sum(1);obs=-np.log1p(-np.minimum(c[k:]**2,1-1e-14)).sum()
  ps.append((1+(stats>=obs).sum())/20001)
 return dict(pairwise_correlations=G,canonical_correlations=c,closed_stepwise_p=np.maximum.accumulate(ps))
primary=df.loc[df.selected_primary,'accession'].tolist()
sets={'primary':primary,'no_source_cautions':df.loc[df.selected_primary & ~df.source_caution_flag,'accession'].tolist(),'plus_model_peptides':df.loc[(df.selected_primary | (df.category=='controlled_model_context')) & ~df.missing_value_warning,'accession'].tolist(),'plus_molecular_descriptors':primary+df.loc[df.category=='constructed_molecular_descriptor','accession'].tolist()}
res={'selection_sha256':hashlib.sha256((P/'selection_manifest.json').read_bytes()).hexdigest(),'seed':20261005,'PA_replicates':B,'CCA_permutations':20000,'models':{}}
for name,codes in sets.items():
 M=centered(np.column_stack([D[c] for c in codes]));Z=M/M.std(0,ddof=1)
 variants={'standardized':Z}
 if name=='primary':
  dist=np.maximum(0,1-np.abs(np.corrcoef(Z.T)));np.fill_diagonal(dist,0)
  groups=fcluster(linkage(squareform(dist,checks=False),method='complete'),.1,criterion='distance')
  W=Z/np.sqrt(np.array([sum(groups==g) for g in groups]));variants.update(raw=M,cluster_balanced=W)
  pd.DataFrame({'accession':codes,'cluster':groups}).to_csv(P/'redundancy_clusters.csv',index=False)
  pd.DataFrame(np.column_stack([D[c] for c in codes]),index=order,columns=codes).rename_axis('aa').to_csv(P/'selected_intrinsic_values.csv')
 for v,x in variants.items():
  key=name+'_'+v;u,s,vt=decomp(x);p=pa(x);item={'n_indices':len(codes),'codes':codes,'spectrum_percent':100*s*s/(s@s),'parallel_analysis':p,'cross_3PC':cross(u[:,:3]),'cross_retained':cross(u[:,:max(1,p['retained'])]),'matrix_rank':int(np.linalg.matrix_rank(x))}
  if v=='cluster_balanced':item['block_parallel_analysis']=pa(x,groups)
  if name=='primary' and v=='standardized':
   item['top_selected_correlations']={str(k+1):sorted([{'code':c,'description':R[next(i for i,r in enumerate(R) if r['accession']==c)]['fields']['D'],'r':float(np.corrcoef(u[:,k],D[c])[0,1])} for c in codes],key=lambda a:abs(a['r']),reverse=True)[:8] for k in range(5)}
   item['leave_one_index_spectrum_range']=np.array([np.min([100*decomp(np.delete(x,j,axis=1))[1][:5]**2/np.sum(decomp(np.delete(x,j,axis=1))[1]**2) for j in range(x.shape[1])],axis=0),np.max([100*decomp(np.delete(x,j,axis=1))[1][:5]**2/np.sum(decomp(np.delete(x,j,axis=1))[1]**2) for j in range(x.shape[1])],axis=0)])
  pd.DataFrame(u[:,:5],index=order,columns=[f'PC{i}' for i in range(1,6)]).rename_axis('aa').to_csv(P/(key+'_coordinates.csv'))
  res['models'][key]=item; print(key,'n',len(codes),'spectrum',np.round(item['spectrum_percent'][:5],2),'retained',p['retained'],'CCA',item['cross_3PC'],flush=True)
  (P/'results.json').write_text(json.dumps(res,default=lambda a:a.tolist() if isinstance(a,np.ndarray) else a.item(),indent=2))
