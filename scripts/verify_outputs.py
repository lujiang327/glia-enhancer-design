#!/usr/bin/env python3
"""Check measured output arithmetic and snapshot artifact checksums."""
import csv
import hashlib
import json
from pathlib import Path


def rows(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))


def run():
    q=Path('reports/qc');a=json.loads((q/'fragment_audit.json').read_text());p=json.loads((q/'pool_qc.json').read_text());d=json.loads((q/'donor_qc.json').read_text())
    assert p['counts']['input']==sum(v for k,v in p['counts'].items() if k!='input')==a['fragment_records']
    assert d['flow']['input']==sum(v for k,v in d['flow'].items() if k!='input')
    ds=rows(q/'donor_summary.tsv');assert sum(int(r['retained_fragment_records']) for r in ds)==d['flow']['retained']
    cells=rows(q/'cells.unambiguous.tsv');assert len(cells)==len({r['cell_id'] for r in cells})==len({r['bare_barcode'] for r in cells})==3434
    f=rows(q/'frip.tsv');byname={r['scope']:r for r in f}
    for r in ds:assert int(byname[r['donor']]['retained_fragment_records'])==int(r['retained_fragment_records'])
    for r in f:
        n,k=int(r['retained_fragment_records']),int(r['fragments_overlapping_peak']);assert 0<=k<=n;assert abs(float(r['FRiP'])-k/n)<1e-12
    assert int(byname['full_pool']['retained_fragment_records'])==p['counts']['retained']
    paths=[]
    for root in ['scripts','config','reports','data/metadata','data/raw','data/intermediate/donor_pseudobulks']:
        paths.extend(p for p in Path(root).rglob('*') if p.is_file() and '__pycache__' not in str(p))
    with open('logs/artifact_sha256.tsv','w') as out:
        w=csv.writer(out,delimiter='\t');w.writerow(['path','bytes','sha256'])
        for path in sorted(paths):
            h=hashlib.sha256()
            with path.open('rb') as f:
                for block in iter(lambda:f.read(2**20),b''):h.update(block)
            w.writerow([str(path),path.stat().st_size,h.hexdigest()])
    print('PASS: fragment accounting, unique mapped cell IDs, donor counts, FRiP arithmetic, and artifact checksums.')


if __name__=='__main__':run()
