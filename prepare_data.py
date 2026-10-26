"""Restore the frozen data beside this source tree from the Zenodo data ZIP.

Usage: python prepare_data.py /path/to/ASAI_Zenodo_data_and_code_2026-10-26.zip
"""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import shutil
from zipfile import ZipFile

DATA_ROOTS={'inputs','work','outputs','source_data','intrinsic','rsa_pc3','kingdom_null','manuscript_tables','selected_figure_sources'}
CODE_SUFFIXES={'.py','.js','.cjs','.r'}
EXPECTED_RSA_SHA256='412039f02ff0fb37730e41eacd1bda6d9defb0922318d927e46088c8a0b147f5'

def restore(archive,root):
    root.mkdir(parents=True,exist_ok=True)
    restored=0
    with ZipFile(archive) as z:
        members=[m for m in z.infolist() if not m.is_dir()]
        marker=[m.filename for m in members if m.filename.endswith('/inputs/residue_RSA_6908_corrected.csv')]
        if len(marker)!=1:raise ValueError('Expected exactly one frozen residue dataset in the archive.')
        prefix=marker[0].split('/inputs/')[0]+'/'
        for member in members:
            if not member.filename.startswith(prefix):continue
            relative=PurePosixPath(member.filename[len(prefix):])
            if relative.is_absolute() or '..' in relative.parts:raise ValueError('Unsafe archive path.')
            if not relative.parts or relative.parts[0] not in DATA_ROOTS:continue
            if relative.suffix.lower() in CODE_SUFFIXES:continue
            target=root.joinpath(*relative.parts)
            if root.resolve() not in target.resolve().parents:raise ValueError('Unsafe extraction destination.')
            target.parent.mkdir(parents=True,exist_ok=True)
            with z.open(member) as source,target.open('wb') as destination:shutil.copyfileobj(source,destination)
            restored+=1
    residue=root/'inputs/residue_RSA_6908_corrected.csv'
    if hashlib.sha256(residue.read_bytes()).hexdigest()!=EXPECTED_RSA_SHA256:
        raise ValueError('The residue input does not match the frozen analysis dataset.')
    print(f'Restored {restored} data and result files; frozen RSA checksum verified.')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args()
    restore(args.archive,args.root.resolve())
