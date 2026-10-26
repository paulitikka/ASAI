"""Export the six tables and their full-precision sources from one calculation."""
from pathlib import Path
import csv,json,math
from decimal import Decimal,InvalidOperation

P=Path(__file__).resolve().parent
AA=list('ACDEFGHIKLMNPQRSTVWY')

def load(path):return json.loads(path.read_text(encoding='utf-8'))
def csv_rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(path,rows,fields=None,bom=False):
    path.parent.mkdir(parents=True,exist_ok=True)
    if fields is None:fields=list(rows[0])
    with path.open('w',encoding='utf-8-sig' if bom else 'utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
def rsa_interval(v):
    if v==0:return 0
    k=round(v*10)
    if 1<=k<=9 and abs(v-k/10)<=1e-12:v=k/10
    return next((i for i in range(1,10) if v<=i/10),10)
def raw_interval(v):
    if v==0:return 0
    return next((i for i in range(1,10) if v<=i/10),10)

def main():
    out=P/'outputs';src=P/'source_data';tables=out/'tables'
    d=load(out/'decompose.json');m=load(out/'admissible_parameters.json');windows=load(out/'s6_windows.json')
    rows=d['table1FullPrecision'];joint=d['table2FullPrecision'];s1=d['S5FullPrecision']
    write(src/'Table1_full_precision.csv',rows)
    write(src/'TableS2_full_precision.csv',joint)
    write(src/'TableS1_full_precision.csv',s1)
    scores=load(out/'model_scores_by_split.json')
    full_models=[]
    for key in m['models']:
        values=[r[key] for r in scores];mean=math.fsum(values)/len(values)
        full_models.append({'model':key,'mean_held_out_R2':mean,'unexplained_percent':100*(1-mean),'n_type_split_evaluations':len(values)})
    write(src/'Table3_full_precision.csv',full_models)
    write(src/'Table3_scores_by_split_type.csv',[{'split':r['split']+1,'aa':r['type'],**{k:r[k] for k in m['models']}} for r in scores])
    write(src/'Table3_fit_diagnostics.csv',load(out/'model_diagnostics.json'))
    full_windows=[{'comparison':r['label'],'N_euk':r['nEuk'],'N_prok':r['nProk'],
                   'mean_length_euk':r['meanLenEukFullPrecision'],'mean_length_prok':r['meanLenProkFullPrecision'],
                   'PC1_t_euk_minus_prok':r['pc1_tFullPrecision'],'PC2_t_euk_minus_prok':r['pc2_tFullPrecision']} for r in windows]
    write(src/'TableS3_full_precision.csv',full_windows)
    write(tables/'Table_1.csv',[{'Type':r['t'],'Mean RSA':f"{r['meanRSA']:.4f}",**{k:f'{r[k]:+.4f}' for k in ['PC1','PC2','PC3']}} for r in rows],bom=True)
    write(tables/'Table_2.csv',[{'#':r['n'],'Parameter':'Signed charge balance' if r['n']==4 else r['name'],'Terms in the model':'1 (log10)' if r['n']==1 else r['terms']} for r in m['parameters']],bom=True)
    names={'length':'Length alone','composition':'Length + composition','charge':'Length + charge + signed charge balance','hydrophobicity':'Length + sequence-encoded hydrophobicity','all':'All five together'}
    write(tables/'Table_3.csv',[{'Model':names[r['model']],'Mean held-out R²':f"{r['mean_held_out_R2']:.4f}",'Unexplained, %':f"{r['unexplained_percent']:.1f}%"} for r in full_models],bom=True)
    write(tables/'Table_S1.csv',[{'RSA interval, %':'>90' if r['interval']=='90+' else r['interval'],
        'Chain-type pairs':r['pairs'],'% of all':f"{r['pct']:.2f}",'CV across types':f"{r['cv']:.3f}"} for r in s1],bom=True)
    columns={'PC1 euk':'PC1_euk','PC1 prok':'PC1_prok','ΔPC1':'dPC1','PC2 euk':'PC2_euk','PC2 prok':'PC2_prok','ΔPC2':'dPC2'}
    write(tables/'Table_S2.csv',[{'Type':r['t'],**{title:f'{r[key]:+.4f}' for title,key in columns.items()}} for r in joint],bom=True)
    write(tables/'Table_S3.csv',[{'Comparison':r['comparison'],'N euk':r['N_euk'],'N prok':r['N_prok'],
        'mean length euk, residues':f"{r['mean_length_euk']:.1f}",'mean length prok, residues':f"{r['mean_length_prok']:.1f}",
        'PC1 t':f"{r['PC1_t_euk_minus_prok']:+.2f}",'PC2 t':f"{r['PC2_t_euk_minus_prok']:+.2f}"} for r in full_windows],bom=True)
    fractions=[]
    for i,interval in enumerate(s1):
        for j,aa in enumerate(d['aminoAcidOrder']):
            count=d['intervalCountsByType'][i][j];den=d['chainsWithType'][j]
            fractions.append({'interval_order':i,'interval':interval['interval'],'aa':aa,'pairs_in_interval':count,'chains_with_type':den,'fraction_within_type':count/den})
    write(src/'TableS1_fractions_by_type.csv',fractions)

    # Export one canonical chain/type aggregation for every downstream analysis.
    pool=load(P/'work/pooled-acc.json');assignment=load(out/'assign.json')
    selected=csv_rows(P/'inputs/selected_6908_organisms.csv')
    labels={}
    for r in selected:
        key=(r['pdb']+r['chain']).lower()
        if key in labels:raise ValueError('Nonunique legacy chain key.')
        labels[key]=(r['pdb'].lower(),r['chain'])
    order='ARNDCQEGHILKMFPSTWYV';chain_rows=[];events=[]
    for key,c in pool.items():
        pdb,chain=labels[key]
        for j,aa in enumerate(order):
            if not c['cnt'][j]:continue
            value=c['sum'][j]/c['cnt'][j]
            chain_rows.append({'pdb':pdb,'chain':chain,'aa':aa,'mean_RSA':value,'n_residues':c['cnt'][j],'chain_key':key})
            old,new=raw_interval(value),rsa_interval(value)
            if old!=new:events.append({'pdb':pdb,'chain':chain,'aa':aa,'unmodified_chain_mean_RSA':repr(value),'old_interval':old,'new_interval':new,'nearest_boundary':round(value*10)/10})
    chain_rows.sort(key=lambda r:(r['pdb'],r['chain'],r['aa']))
    write(out/'chain_amino_acid_mean_RSA.csv',chain_rows)
    write(out/'interval_boundary_adjustments.csv',events,fields=['pdb','chain','aa','unmodified_chain_mean_RSA','old_interval','new_interval','nearest_boundary'])
    annotations=[]
    for r in selected:
        key=(r['pdb']+r['chain']).lower()
        annotations.append({'pdb':r['pdb'],'chain':r['chain'],'sequence':r['sequence'],'organism':r['organism'],
            'legacy_taxid_text':r['taxid'],'entity_id':r['entity_id'],'kingdom_classification':assignment[key]['k'],
            'n_RSA_residues':pool.get(key,{}).get('L',0),'has_processed_RSA':key in pool})
    write(out/'chain_selection_organism_classification.csv',annotations)
    (src/'pca_variance.json').write_text(json.dumps({'all_types_percent':d['shares1FullPrecision'],'joint_kingdom_percent':d['shares2FullPrecision']},indent=2)+'\n')
    # Printed Word tables are reference-only and are opened after all exports.
    dif=[]
    for name in ['1','2','3','S1','S2','S3']:
        ref=P/'reference/printed_Final13'/f'Table_{name}.csv'
        if not ref.exists():continue
        old=csv_rows(ref);new=csv_rows(tables/f'Table_{name}.csv')
        for a,b in zip(old,new):
            for (title,x),y in zip(a.items(),b.values()):
                try:nx=Decimal(x.replace('−','-').replace('%','').strip());ny=Decimal(str(y).replace('−','-').replace('%','').strip())
                except InvalidOperation:continue
                if nx!=ny:dif.append({'table':name,'row':next(iter(a.values())),'column':title,'Word':x,'r3':y})
    write(out/'differences_from_printed_Final13.csv',dif)
    (out/'export_summary.json').write_text(json.dumps({'n_tables':6,'n_chain_type_pairs':len(chain_rows),
        'n_interval_boundary_reassignments':len(events),'n_numeric_differences_from_printed_Final13':len(dif),
        'rounding':'Full-precision calculation; one final rounding at display export; no double rounding',
        'interval_rule':'Exact zeros separate; (lower,upper]; values within 1e-12 of a positive 0.1 boundary treated as that boundary only for interval assignment'},indent=2)+'\n')
    print('Exported six tables, full-precision source data and boundary audit.')

if __name__=='__main__':main()
