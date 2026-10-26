from pathlib import Path
import json,hashlib,csv
P=Path(__file__).resolve().parent
f=P/'inputs/residue_RSA_6908_corrected.csv'
assert hashlib.sha256(f.read_bytes()).hexdigest()=='412039f02ff0fb37730e41eacd1bda6d9defb0922318d927e46088c8a0b147f5'
ids=list(csv.DictReader(open(P/'inputs/intrinsic_51_manifest.csv')));assert len(ids)==51 and len({x['accession'] for x in ids})==51
D=json.load(open(P/'outputs/decompose.json'));assert D['totalPairs']==120822 and D['kingdomCounts']=={'euk':1778,'prok':2877}
M=json.load(open(P/'outputs/admissible_parameters.json'));assert M['chains']==6287 and M['models']['all']['meanHeldOutR2']==.1358
N=json.load(open(P/'kingdom_null/results.json'));assert N['half_split']['all_same_sign_count']==348 and N['original_sizes']['mean_shift_extreme_count']==0
V=json.load(open(P/'rsa_pc3/validation_results.json'));assert V['best_intrinsic']['accession']=='FAUJ880112' and abs(V['best_intrinsic']['p_max_51']-.8545072746362682)<1e-12
print('PASS: input digest, 51 unique indices, chain groups/pairs, full model, protein null and PC3 multiplicity.')
