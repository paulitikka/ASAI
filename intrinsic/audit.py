from pathlib import Path
import json,re,hashlib
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import squareform
P=Path(__file__).resolve().parent
R=json.loads((P.parent/'inputs/aaindex_records.json').read_text()); order=list('ARNDCQEGHILKMFPSTWYV')
# Explicit record decisions frozen before inspecting new spectra. Positions refer to bundled JSON.
core=[13,14,20,21,22,*range(25,31),58,*range(71,77),*range(78,90),95,109,110,111,113,115,116,117,132,144,145,151,176,239,*range(313,318),354,355,382,387,388,396,398,399,400,401,506]
oper=[11,12,15,16,17,32,177,178,179,180,238,433,477,478,479,480,482,512,513]
graph=[77,378,475,*range(531,553)]
conf=[90,91,92,93,94,255,256,356,362,363,375,376,377,422,423]
unc=[0,8,62,68,69,70,107,108,131,149,*range(152,159),251,318,320,351,357,380,397,514,529,530]
composite=[67,*range(358,362),389,390,391,508]
assert len(set(core+oper+graph+conf+unc+composite+[211]))==len(core+oper+graph+conf+unc+composite+[211])
rows=[]
for i,r in enumerate(R):
 f=r['fields'];j=f.get('J','');cat='excluded_context';reason='Protein/sequence/structure statistics, fitted structural propensity, or biological-context endpoint; see description and source title.'
 if i in core:cat='intrinsic_candidate';reason='Chemical composition, geometry, charge, or physicochemical measurement/calculation independent of native-protein statistics; judged from record and cited provenance.'
 if i in oper:cat='controlled_model_context';reason='Controlled peptide measurement or model-tripeptide geometry; operationally context-dependent, included only in expanded sensitivity.'
 if i in graph:cat='constructed_molecular_descriptor';reason='Protein-independent molecular graph/EIIP descriptor; legitimate molecular information but representation-dependent; separate sensitivity.'
 if i in conf:cat='conformational_model';reason='Helix/coil or conformer energetics/populations depend on a polymer or conformational model; excluded from primary chemistry set.'
 if i in unc:cat='uncertain_provenance';reason='Potential chemical property, but independent origin of these exact values is insufficiently established; excluded rather than presumed intrinsic.'
 if i in composite:cat='excluded_composite';reason='Previously combined, optimized, or principal-component descriptor; not a primitive measurement for estimating dimensionality.'
 if i==211:cat='excluded_quality';reason='Record explicitly reports missing values filled with zeros; do not interpret placeholders as measurements.'
 missing=bool(re.search('missing|filled.*zero',j,re.I));flag=bool(re.search(r'!|missing|filled.*zero|\((?:Pro|Arg)\s+[0-9]',j,re.I))
 selected=i in core and not missing
 rows.append(dict(accession=r['accession'],description=f.get('D',''),category=cat,selected_primary=selected,missing_value_warning=missing,source_caution_flag=flag,reason=reason,source_title=f.get('T',''),source_citation=j,source_reference=f.get('R',''),evidence='AAindex record and provenance judgment; selected targeted source abstracts checked; not a full original-table numeric audit'))
df=pd.DataFrame(rows);df.to_csv(P/'AAindex_553_audit.csv',index=False)
manifest={'primary':df.loc[df.selected_primary,'accession'].tolist(),'definition':'Protein-independent chemistry; exclude explicit missing-value placeholders. No guarantee of context-free universal values.','categories':df.category.value_counts().to_dict()}
(P/'selection_manifest.json').write_text(json.dumps(manifest,indent=2));print('FROZEN',hashlib.sha256((P/'selection_manifest.json').read_bytes()).hexdigest(),len(manifest['primary']),flush=True)
