"""Regenerate downstream RSA coordinates from the current canonical chain means."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd

P=Path(__file__).resolve().parent
order=list('ARNDCQEGHILKMFPSTWYV')
c=pd.read_csv(P/'outputs/chain_amino_acid_mean_RSA.csv',keep_default_na=False,float_precision='round_trip')
a=c.pivot(index=['pdb','chain'],columns='aa',values='mean_RSA').reindex(columns=order)
q=np.nanquantile(a.to_numpy(),np.arange(1,100)/100,axis=0).T
u,s,v=np.linalg.svd(q-q.mean(0),full_matrices=False)
ref=pd.DataFrame(json.loads((P/'outputs/decompose.json').read_text())['table1FullPrecision']).set_index('t').loc[order]
for k in range(6):
    if k<3:sign=np.dot(u[:,k],ref[f'PC{k+1}'])
    else:sign=u[np.argmax(abs(u[:,k])),k]
    if sign<0:u[:,k]*=-1
if np.max(np.abs(u[:,:3]-ref[['PC1','PC2','PC3']].to_numpy()))>5e-9:raise ValueError('Independent SVD disagrees with canonical Jacobi PCA.')
frame=pd.DataFrame({'aa':order,'mean_chain_RSA':a.mean().to_numpy(),'sd_chain_RSA':a.std(ddof=1).to_numpy(),
                    **{f'RSA_PC{k+1}':u[:,k] for k in range(6)}})
for folder in ['intrinsic','rsa_pc3']:
    frame.to_csv(P/folder/'rsa_coordinates.csv',index=False)
hashes={'source':'Current outputs/chain_amino_acid_mean_RSA.csv; independent SVD matches the JS PCA',
        'residue_sha256':hashlib.sha256((P/'inputs/residue_RSA_6908_corrected.csv').read_bytes()).hexdigest(),
        'canonical_chain_means_sha256':hashlib.sha256((P/'outputs/chain_amino_acid_mean_RSA.csv').read_bytes()).hexdigest(),
        'rsa_coordinates_sha256':hashlib.sha256((P/'intrinsic/rsa_coordinates.csv').read_bytes()).hexdigest()}
(P/'outputs/downstream_coordinate_provenance.json').write_text(json.dumps(hashes,indent=2)+'\n')
print('Canonical downstream RSA coordinates regenerated and checked independently.')

# Rebuild the frozen selection audit before the PC3 multiple-testing scan.
import subprocess,sys
subprocess.run([sys.executable,str(P/'intrinsic/audit.py')],check=True,cwd=P)
