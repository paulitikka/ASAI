"""Recalculate Tables 1-3 and S1-S3 using Python's standard library and Node."""
import sys
from runtime import ROOT,verify_inputs,locate_node,step,environment

def main():
    environment();verify_inputs()
    with (ROOT/'outputs/table_run.log').open('w',encoding='utf-8') as log:
        step('1/3 Calculating the six tables',[locate_node(),'run_js.cjs','--tables'],log)
        step('2/3 Exporting tables and full-precision source data',[sys.executable,'export_results.py'],log)
        step('3/3 Verifying table consistency',[sys.executable,'verify_tables.py'],log)
    print('PASS. Tables: outputs/tables/. Verification: outputs/table_verification.json',flush=True)

if __name__=='__main__':main()
