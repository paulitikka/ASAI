"""Verify the six tables, precision, model diagnostics and interval arithmetic."""
from pathlib import Path
import csv,json,math

P=Path(__file__).resolve().parent
def load(p):return json.loads(p.read_text())
def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def close(a,b,tol=1e-11):
    if abs(float(a)-float(b))>tol:raise AssertionError((a,b))
def mean(a):return math.fsum(a)/len(a)
def sd(a):
    m=mean(a);return math.sqrt(math.fsum((v-m)**2 for v in a)/(len(a)-1))

def main():
    d=load(P/'outputs/decompose.json');model=load(P/'outputs/admissible_parameters.json');w=load(P/'outputs/s6_windows.json')
    assert d['totalPairs']==120822 and d['kingdomCounts']=={'euk':1778,'prok':2877}
    assert model['chains']==6287
    t1=d['table1FullPrecision'];t2=d['table2FullPrecision']
    for pc in ['PC1','PC2','PC3']:close(math.fsum(r[pc]**2 for r in t1),1)
    close(math.fsum(r['PC1']*r['PC2'] for r in t1),0)
    for pc in ['PC1','PC2']:
        close(math.fsum(r[pc+'_euk']**2+r[pc+'_prok']**2 for r in t2),1)
        for r in t2:close(r['d'+pc],r[pc+'_euk']-r[pc+'_prok'])
    assert len(d['intervalCountsByType'])==11
    for j,total in enumerate(d['chainsWithType']):assert sum(r[j] for r in d['intervalCountsByType'])==total
    for i,r in enumerate(d['S5FullPrecision']):
        fractions=[c/n for c,n in zip(d['intervalCountsByType'][i],d['chainsWithType'])]
        assert sum(d['intervalCountsByType'][i])==r['pairs']
        close(r['pct'],100*r['pairs']/120822)
        close(r['cv'],sd(fractions)/mean(fractions))
    assert sum(r['pairs'] for r in d['S5FullPrecision'])==120822
    close(sum(r['pct'] for r in d['S5FullPrecision']),100)
    scores=load(P/'outputs/model_scores_by_split.json');diag=load(P/'outputs/model_diagnostics.json')
    assert len(scores)==400 and len(diag)==2000
    assert len({(r['split'],r['type']) for r in scores})==400
    for r in diag:close(r['R2'],1-r['SSE']/r['SST'])
    for key,m in model['models'].items():
        values=[r[key] for r in scores];assert m['n']==len(values)==400
        assert f'{mean(values):.4f}'==f"{m['meanHeldOutR2']:.4f}"
        assert f'{100*(1-mean(values)):.1f}'==f"{m['unexplainedPct']:.1f}"
    for r in w:
        for pc in ['PC1','PC2']:
            values=r['delta'+pc+'FullPrecision'];assert len(values)==20
            close(r[pc.lower()+'_tFullPrecision'],mean(values)/(sd(values)/math.sqrt(20)))
    expected_sizes={'1':20,'2':5,'3':5,'S1':11,'S2':20,'S3':3}
    for name,n in expected_sizes.items():assert len(rows(P/'outputs/tables'/f'Table_{name}.csv'))==n
    printed=rows(P/'outputs/tables/Table_S1.csv')
    for actual,r in zip(printed,d['S5FullPrecision']):
        assert int(actual['Chain-type pairs'])==r['pairs']
        assert actual['% of all']==f"{r['pct']:.2f}"
        assert actual['CV across types']==f"{r['cv']:.3f}"
    for actual,r in zip(rows(P/'outputs/tables/Table_1.csv'),t1):
        assert actual['Mean RSA']==f"{r['meanRSA']:.4f}"
        for key in ['PC1','PC2','PC3']:assert actual[key]==f'{r[key]:+.4f}'
    check={'status':'PASS','version':'2026.10.26-r3','n_tables':6,'n_chain_type_pairs':120822,
           'n_model_type_split_scores':400,'n_model_fit_diagnostics':2000,
           'checks':['PCA unit vectors and directions','joint coordinate differences',
                     'type-specific interval denominators and CV','percentages and single final rounding',
                     'held-out SSE/SST and model averaging','full-precision signed paired t statistics'],
           'scope':'Arithmetic and export consistency; full independent raw-data checks are in verify_independent.py.'}
    (P/'outputs/table_verification.json').write_text(json.dumps(check,indent=2)+'\n')
    print('PASS: all six tables, full-precision source values and model diagnostics are consistent.')

if __name__=='__main__':main()
