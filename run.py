"""Reproduce the processed-data analysis. Use --full for every article analysis."""
import argparse,json,sys
from runtime import ROOT,verify_inputs,locate_node,step,environment

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full',action='store_true',help='Also rerun 200,000 whole-chain allocations and intrinsic-AA PA/CCA.')
    parser.add_argument('--no-figures',action='store_true')
    args=parser.parse_args();environment();verify_inputs()
    stages=[('Calculating RSA, context models and corrected permutation tests',[locate_node(),'run_js.cjs']),
            ('Exporting six tables and numerical source data',[sys.executable,'export_results.py']),
            ('Preparing canonical RSA coordinates for downstream analyses',[sys.executable,'prepare_analysis.py']),
            ('Rebuilding PC3 distribution descriptors',[sys.executable,'rsa_pc3/explore.py']),
            ('Validating PC3 asymmetry and property correlations',[sys.executable,'rsa_pc3/validate.py'])]
    if args.full:
        stages += [('Running 200,000 whole-chain random allocations',[sys.executable,'kingdom_null/test.py']),
                   ('Running intrinsic-AA parallel analysis and CCA',[sys.executable,'intrinsic/analyze.py']),
                   ('Running intrinsic-AA family sensitivity',[sys.executable,'intrinsic/sensitivity.py']),
                   ('Exporting intrinsic PC2 interpretation',[sys.executable,'intrinsic/interpret_pc2.py'])]
    if not args.no_figures:stages.append(('Drawing Figures 1-4',[sys.executable,'plot_figures.py']))
    stages.append(('Verifying the six tables',[sys.executable,'verify_tables.py']))
    if args.full:stages.append(('Independently checking raw data, PCA and all 2000 model fits',[sys.executable,'verify_independent.py']))
    with (ROOT/'outputs/run.log').open('w',encoding='utf-8') as log:
        for i,(title,command) in enumerate(stages,1):step(f'{i}/{len(stages)} {title}',command,log)
    (ROOT/'outputs/run_scope.json').write_text(json.dumps({'version':'2026.10.26-r3','full':args.full,'figures':not args.no_figures,'all_stages_completed':True},indent=2)+'\n')
    print('PASS. Analysis completed. Results: outputs/, source_data/, rsa_pc3/, intrinsic/, kingdom_null/.',flush=True)

if __name__=='__main__':main()
