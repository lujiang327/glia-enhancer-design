#!/usr/bin/env python3
"""Verify the UCSC analysis-set download and create a canonical hg38 FASTA."""
import gzip
import hashlib
import json
from pathlib import Path
from pyfaidx import Fasta


def digest(path, algorithm='sha256'):
    h=hashlib.new(algorithm)
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()


def run():
    root=Path('data/raw/reference');source=root/'hg38.analysisSet.fa.gz'
    if not source.exists():source=source.with_suffix(source.suffix+'.part')
    expected=dict((name,value) for value,name in (l.split() for l in (root/'upstream.md5sum.txt').read_text().splitlines()))['hg38.analysisSet.fa.gz']
    assert digest(source,'md5')==expected,'Reference download MD5 mismatch'
    if source.name.endswith('.part'):source=source.rename(root/'hg38.analysisSet.fa.gz')
    canonical={f'chr{i}' for i in range(1,23)}|{'chrX'}
    target=root/'hg38.canonical.fa';partial=target.with_suffix('.fa.part');lengths={};name=None
    with gzip.open(source,'rt') as f,partial.open('w') as out:
        for line in f:
            if line.startswith('>'):
                name=line[1:].split()[0]
                if name in canonical:
                    assert name not in lengths
                    lengths[name]=0;out.write('>'+name+'\n')
            elif name in canonical:
                seq=line.strip().upper();assert set(seq)<=set('ACGTN')
                lengths[name]+=len(seq);out.write(seq+'\n')
    upstream=dict((c,int(n)) for c,n in (l.split() for l in Path('data/raw/hg38.chrom.sizes').read_text().splitlines()))
    assert set(lengths)==canonical
    assert all(lengths[c]==upstream[c] for c in canonical),'Reference length mismatch'
    partial.replace(target)
    with Fasta(str(target),rebuild=True):pass
    order=sorted(canonical)
    (root/'hg38.canonical.chrom.sizes').write_text(''.join(f'{c}\t{lengths[c]}\n' for c in order))
    report=dict(source_url='https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz',source_md5=expected,source_sha256=digest(source),fasta=str(target),fasta_sha256=digest(target),chromosomes=lengths,transformation='Uppercase chr1-22,X only; coordinates and hard-masked N bases preserved.')
    out=Path('reports/preparation');out.mkdir(exist_ok=True)
    (out/'reference.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':run()
