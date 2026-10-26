"""Reproduce the descriptive PC2 interpretation; no new dimension-selection test."""
from pathlib import Path
import json,html
import numpy as np
import pandas as pd
P=Path(__file__).resolve().parent
x=pd.read_csv(P/'selected_intrinsic_values.csv',index_col=0)
z=(x-x.mean())/x.std(ddof=1)
u,s,vt=np.linalg.svd(z.values,full_matrices=False)
# Choose signs for readable chemical interpretation; PCA signs are arbitrary.
if np.corrcoef(u[:,0],x['FAUJ830101'])[0,1]<0:u[:,0]*=-1
if np.corrcoef(u[:,1],x['FASG760101'])[0,1]<0:u[:,1]*=-1
corr=lambda a,b:float(np.corrcoef(a,b)[0,1])
props={'FAUJ880104':'Side-chain length (STERIMOL)','CHAM830106':'Bonds in longest side-chain path','FASG760101':'Molecular weight','FAUJ880103':'Side-chain van der Waals volume','FAUJ880109':'Hydrogen-bond donor count','WOLR810101':'Vapor-to-water hydration free energy'}
rows=[];design=np.column_stack([np.ones(20),u[:,0]])
for code,label in props.items():
 y=z[code].values;res=y-design@np.linalg.lstsq(design,y,rcond=None)[0]
 rows.append({'accession':code,'property':label,'r_PC2':corr(y,u[:,1]),'r_PC2_after_removing_PC1':corr(res,u[:,1])})
t=pd.DataFrame(rows);t.to_csv(P/'PC2_supporting_correlations.csv',index=False)
X=np.column_stack([np.ones(20),z[['FASG760101','FAUJ830101']].values]);fit=X@np.linalg.lstsq(X,u[:,1],rcond=None)[0]
r2=1-float(np.sum((u[:,1]-fit)**2)/np.sum(u[:,1]**2))
result={'n_AA':20,'n_indices':51,'PC1_PC2_correlation':corr(u[:,0],u[:,1]),'variance_percent':(100*s[:2]**2/np.sum(s*s)).tolist(),'correlations':rows,'PC2_mass_plus_hydrophobicity_R2':r2,'interpretation':'PC1: Hydrophobicity; PC2: Side-chain size. Descriptive labels, not pure physical quantities or causal independence.','status':'Post hoc interpretation of fixed primary PCA; not independent validation. No rotation performed.'}
(P/'PC2_interpretation_results.json').write_text(json.dumps(result,indent=2))
paragraphs=[
'The agreed descriptive labels are PC1: Hydrophobicity and PC2: Side-chain size. PC1 and PC2 explain 32.07% and 27.34% of the standardized 51-index matrix, respectively.',
'Proposed manuscript wording: “The second principal component primarily represents side-chain size variation orthogonal to PC1, with additional contributions from hydration and hydrogen bonding.”',
'PC2 correlates with side-chain length (r = 0.854), the number of bonds in the longest side-chain path (r = 0.851), molecular weight (r = 0.794), and van der Waals volume (r = 0.729). These consistent associations support its interpretation as a size/extension axis.',
'After subtracting the linear component associated with PC1 from each property, residual van der Waals volume correlates with PC2 at r = 0.952 and residual molecular weight at r = 0.926. Because PC2 is orthogonal to PC1, these are also the partial correlations controlling for PC1. “At fixed hydrophobicity” here means controlling for the hydrophobicity-associated PC1, not fixing every experimental hydrophobicity scale.',
f'A linear model containing molecular weight and the Fauchere–Pliska experimental hydrophobicity index reconstructs {100*r2:.1f}% of PC2 variation across the 20 amino acids. This is descriptive reconstruction: both predictors contributed to the PCA, so it is not independent predictive validation.',
'PC2 also correlates with hydrogen-bond donor count (r = 0.815) and vapor-to-water hydration free energy (r = −0.865). Arginine, lysine and histidine have high scores, while glycine and alanine have low scores. This explains why “size” is a useful shorthand rather than a complete physical definition of PC2. The sign of a PCA axis is arbitrary; positive PC2 is oriented toward larger molecular weight here.',
'Orthogonality means the two PCA score vectors have zero correlation across these 20 amino acids. It does not require unrelated chemical interpretations: an experimental property can correlate with both axes. Accordingly, the labels avoid calling both axes hydrophobicity, while retaining the hydration contribution in the explanatory text. Neither axis is claimed to be a pure property, and orthogonality does not establish statistical or causal independence.',
'The similar variance shares do not mathematically require a uniquely identifiable physical role for each PC. This interpretation concerns the original, unrotated primary PCA. It does not add a rotation or robustness result beyond the analyses already documented in the bundle.',
'Primary-source context: Wolfenden et al. (1981), “Affinities of amino acid side chains for solvent water,” Biochemistry 20, 849–855, DOI 10.1021/bi00507a030, describes vapor-to-water hydration potentials corrected to pH 7. This thermodynamic comparison differs from partition between water and an organic solvent. Source: https://pubmed.ncbi.nlm.nih.gov/7213619/.',
'Reproduction: run python interpret_pc2.py using the included selected_intrinsic_values.csv. This writes the supporting correlations, numeric results, and this narrative. The primary selection, PCA results and RSA comparisons are unchanged.'
]
(P/'PC2_interpretation.md').write_text('# PC2 interpretation: side-chain size\n\nUpdated 5 October 2026.\n\n'+'\n\n'.join(paragraphs)+'\n\n'+t.round(3).to_csv(index=False))
(P/'PC2_interpretation_section.html').write_text('<h2>PC2 interpretation: side-chain size</h2><p>Added 5 October 2026.</p>'+''.join('<p>'+html.escape(q)+'</p>' for q in paragraphs)+t.round(3).to_html(index=False))
print('PC2 correlations and reconstruction verified:',r2)
