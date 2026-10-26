"""Run from any working directory: python /path/to/run.py [--full]."""
from pathlib import Path
import subprocess,sys,os,json,hashlib
P=Path(__file__).resolve().parent
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
def run(args):subprocess.run(args,cwd=P,check=True)
run([os.environ.get('ASAI_NODE', 'node'),'run_js.cjs'])
# Scientific PC3 follow-up: original algorithms, seed and reference orientation.
import shutil
(P/'rsa_pc3/input').mkdir(parents=True, exist_ok=True)
shutil.copy(P/'inputs/residue_RSA_6908_corrected.csv',P/'rsa_pc3/input/residue_RSA_6908_corrected.csv')
run([sys.executable,'rsa_pc3/explore.py']);run([sys.executable,'rsa_pc3/validate.py'])
if '--full' in sys.argv:
 run([sys.executable,'kingdom_null/test.py'])
 run([sys.executable,'intrinsic/analyze.py'])
run([sys.executable,'export_tables_figures.py'])
print('Completed. New results: outputs/, rsa_pc3/, and (with --full) kingdom_null/, intrinsic/.')
