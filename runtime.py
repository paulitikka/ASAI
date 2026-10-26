from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,time

ROOT=Path(__file__).resolve().parent

def verify_inputs():
    release=json.loads((ROOT/'release.json').read_text())
    for relative,expected in release['inputs'].items():
        p=ROOT/relative
        if not p.is_file():raise SystemExit(f'Missing {relative}; run prepare_data.py with the matching Zenodo data ZIP.')
        if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:raise SystemExit(f'Input checksum mismatch: {relative}')
    design=ROOT/'design/train_test_splits.csv'
    if not design.is_file():raise SystemExit('Missing frozen split design; restore the matching data archive.')
    if hashlib.sha256(design.read_bytes()).hexdigest()!=release['split_design_sha256']:
        raise SystemExit('Split design checksum mismatch.')
    return release

def locate_node():
    if os.environ.get('ASAI_NODE'):return os.environ['ASAI_NODE']
    p=shutil.which('node')
    if p:return p
    if os.name=='nt':
        for var in ['ProgramFiles','ProgramFiles(x86)','LOCALAPPDATA']:
            root=os.environ.get(var)
            if root and (Path(root)/'nodejs/node.exe').is_file():return str(Path(root)/'nodejs/node.exe')
    raise SystemExit('Node.js was not found. Install Node.js 24 and reopen the terminal, or set ASAI_NODE to node.exe.')

def step(title,command,log):
    print(title,flush=True);start=time.monotonic()
    log.write(title+'\n');log.flush()
    result=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:raise SystemExit(f'Calculation stopped at {title}. See {log.name}.')
    log.write(f'Completed in {time.monotonic()-start:.2f} seconds\n');log.flush()

def environment():
    os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
    (ROOT/'outputs').mkdir(exist_ok=True)

