"""Check dataset identity, source-data integrity, figures and an optional manifest."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import numpy as np
import pandas as pd
from pypdf import PdfReader

P=Path(__file__).resolve().parent
AA=set('ACDEFGHIKLMNPQRSTVWY')

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true')
    parser.add_argument('--check-manifest',action='store_true',
        help='Check distributed file bytes before rerunning analyses; regenerated PDF metadata can differ.')
    args=parser.parse_args()
    expected='412039f02ff0fb37730e41eacd1bda6d9defb0922318d927e46088c8a0b147f5'
    assert digest(P/'inputs/residue_RSA_6908_corrected.csv')==expected
    h=pd.read_csv(P/'source_data/Figure1_histogram_bins.csv')
    assert set(h.aa)==AA and len(h)==1200
    assert int(h['count'].sum())==1500829
    integrals=(h.density*(h.bin_right-h.bin_left)).groupby(h.aa).sum()
    assert np.allclose(integrals,1,atol=1e-12)
    f2=pd.read_csv(P/'source_data/Figure2_PC2_dispersion.csv')
    assert len(f2)==20 and set(f2.aa)==AA
    assert -0.9871<float(f2[['sd_chain_mean_RSA','PC2']].corr().iloc[0,1])<-0.986
    f3=pd.read_csv(P/'source_data/Figure3_divergence.csv')
    assert int(f3.n_chain_type_pairs.sum())==120822
    assert np.isclose(f3.percentage_all_pairs.sum(),100)
    assert f3.iloc[-1].interval_label=='>90' and int(f3.iloc[-1].n_chain_type_pairs)==50
    fractions=pd.read_csv(P/'source_data/Figure3_interval_fractions_by_type.csv')
    assert len(fractions)==220 and np.allclose(fractions.groupby('aa').fraction_within_type.sum(),1)
    archived=json.loads((P/'outputs/decompose.json').read_text())
    assert list(f3.n_chain_type_pairs)==[row['pairs'] for row in archived['S5']]
    assert np.allclose(f3.CV_across_types,[row['cv'] for row in archived['S5']],atol=.000501)
    f4=pd.read_csv(P/'source_data/Figure4_kingdom_coordinates.csv')
    assert len(f4)==20 and set(f4.aa)==AA
    assert (f4.delta_PC1_euk_minus_prok>0).sum()==20
    assert (f4.delta_PC2_euk_minus_prok<0).sum()==19
    annotations=pd.read_csv(P/'outputs/chain_selection_organism_classification.csv')
    assert len(annotations)==6908 and int(annotations.has_processed_RSA.sum())==6342
    assert int(annotations.n_RSA_residues.sum())==1500829
    available=annotations[annotations.has_processed_RSA]
    assert (available.kingdom_classification=='euk').sum()==1778
    assert (available.kingdom_classification=='prok').sum()==2877
    stems=['Figure1_residue_RSA','Figure2_PC2_dispersion','Figure3_divergence','Figure4_kingdom_coordinates']
    for stem in stems:
        for ext in ['.pdf','.svg','.png']:
            p=P/'outputs/figures'/(stem+ext)
            assert p.is_file() and p.stat().st_size>1000
        assert len(PdfReader(P/'outputs/figures'/(stem+'.pdf')).pages)==1
    assert len(PdfReader(P/'outputs/figures/Figures_1-4.pdf').pages)==4
    manifest_count=0
    if args.check_manifest:
        assert (P/'MANIFEST.csv').is_file(),'MANIFEST.csv is missing.'
        for row in csv.DictReader((P/'MANIFEST.csv').open()):
            target=P/row['path']
            assert target.is_file(),row['path']
            assert target.stat().st_size==int(row['size_bytes']),row['path']
            assert digest(target)==row['sha256'],row['path']
            manifest_count+=1
    report={'status':'PASS','release_date':'2026-10-23','n_residues':1500829,'n_RSA_chains':6342,
        'n_selected_chains':6908,'n_chain_type_pairs':120822,'n_eukaryotic_chains':1778,
        'n_prokaryotic_chains':2877,'residue_data_sha256':expected,
        'checks':['frozen input checksum','histogram total counts and normalisation',
        'PC2 correlation and twenty type labels','divergence counts and original-pooling consistency',
        'within-type interval fractions','joint kingdom directions','annotation/RSA mapping',
        'four PDF/SVG/PNG figures and four-page collection'],
        'R_renderer':'Not executed; R was unavailable in the packaging environment.'}
    if args.write_report:
        (P/'PACKAGE_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if manifest_count:print(f'PASS: {manifest_count} manifest files match their sizes and SHA256 checksums.')

if __name__=='__main__':main()
