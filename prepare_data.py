"""Restore the matching data deposit without copying saved results into live outputs."""
from pathlib import Path, PurePosixPath
import argparse, hashlib, json, shutil, tempfile, zipfile
from runtime import ROOT, verify_inputs

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path,help='ASAI_Zenodo_data_v2026.10.26-r3.zip')
    args=parser.parse_args()
    version=json.loads((ROOT/'release.json').read_text())['version']
    prefix=f'ASAI_data_v{version}/'
    with zipfile.ZipFile(args.archive) as z:
        members={}
        for item in z.infolist():
            if item.is_dir():continue
            if not item.filename.startswith(prefix):raise SystemExit('Unexpected data archive root.')
            relative=item.filename[len(prefix):]
            path=PurePosixPath(relative)
            if path.is_absolute() or '..' in path.parts or '\\' in relative or relative in members:
                raise SystemExit('Invalid or duplicate archive path.')
            members[relative]=item
        metadata=json.loads(z.read(members['data_release.json']))
        if metadata['version']!=version:raise SystemExit('Data/code release versions differ.')
        manifest=json.loads(z.read(members['SHA256SUMS.json']))
        if set(manifest)!=set(members)-{'SHA256SUMS.json'}:raise SystemExit('Incomplete archive manifest.')
        with tempfile.TemporaryDirectory(prefix='asai-restore-') as tmp:
            stage=Path(tmp)
            for relative,expected in manifest.items():
                payload=z.read(members[relative])
                if hashlib.sha256(payload).hexdigest()!=expected:raise SystemExit(f'Checksum mismatch: {relative}')
                target=stage/relative;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(payload)
            # Validate the entire archive before updating any project file.
            destinations=[]
            for relative in manifest:
                first=PurePosixPath(relative).parts[0]
                if first in {'inputs','design','reference','provenance'}:target=ROOT/relative
                elif first=='numerical_results':target=ROOT/'reference/published_r3'/PurePosixPath(relative).relative_to('numerical_results')
                else:continue
                if target.exists() and target.read_bytes()!=(stage/relative).read_bytes():
                    raise SystemExit(f'Existing file differs; restore into a clean code folder: {target.relative_to(ROOT)}')
                destinations.append((stage/relative,target))
            for source,target in destinations:
                target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
    verify_inputs()
    print(f'Restored and verified data release v{version}.')
    print('Saved numerical results are in reference/published_r3/. New runs write separate live outputs.')

if __name__=='__main__':main()
