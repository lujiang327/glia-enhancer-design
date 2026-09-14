#!/usr/bin/env python3
"""Build molecule-weighted insertion tracks and peak-calling input on a Mac.

Uses validated empirical shifts and retains the existing donor BEDPE files.
Sparse count arrays are retained so pooled and donor BigWigs are reproducible.
"""
import csv
import gzip
import itertools
import json
from pathlib import Path
import numpy as np
import pyBigWig


def write_bigwig(path,sizes,sparse_paths):
    part=path.with_suffix(path.suffix+'.part');bw=pyBigWig.open(str(part),'w');bw.addHeader(list(sizes.items()))
    total=0
    try:
        for c in sizes:
            if len(sparse_paths)==1:
                data=np.load(sparse_paths[0]/f'{c}.npz');pos=data['positions'];values=data['counts']
            else:
                dense=np.zeros(sizes[c],dtype=np.uint32)
                for source in sparse_paths:
                    data=np.load(source/f'{c}.npz');dense[data['positions']]+=data['counts']
                pos=np.flatnonzero(dense).astype(np.int32);values=dense[pos];del dense
            total+=int(values.sum(dtype=np.uint64))
            for k in range(0,len(pos),500000):
                p=pos[k:k+500000];v=values[k:k+500000].astype(float)
                bw.addEntries([c]*len(p),p,ends=p+1,values=v)
    finally:bw.close()
    part.replace(path)
    with pyBigWig.open(str(path)) as check:
        assert int(round(check.header()['sumData']))==total
    return total


def run():
    shift=json.loads(Path('reports/preparation/fragment_shift.json').read_text());assert shift['status']=='PASS'
    delta_a,delta_b=shift['plus_endpoint_delta'],shift['minus_endpoint_delta']
    sizes=dict((c,int(n)) for c,n in (l.split() for l in Path('data/raw/reference/hg38.canonical.chrom.sizes').read_text().splitlines()))
    out=Path('data/intermediate/training');out.mkdir(exist_ok=True);sparse=out/'sparse_counts';sparse.mkdir(exist_ok=True)
    bed=out/'muller.insertions.bed.gz';fragments=out/'muller.original_fragments.tsv.gz';records={};sources=[]
    with open('reports/qc/donor_summary.tsv') as donor_table, gzip.open(str(bed)+'.part','wt',compresslevel=1) as peakinput, gzip.open(str(fragments)+'.part','wt',compresslevel=1) as merged:
        for row in csv.DictReader(donor_table,delimiter='\t'):
            d=row['donor'];dest=sparse/d;dest.mkdir(exist_ok=True);sources.append(dest);seen=set();n=0;endpoint_count=0
            with open(f'data/intermediate/donor_pseudobulks/{d}.bedpe') as f:
                for c,group in itertools.groupby(f,key=lambda l:l.split('\t',1)[0]):
                    assert c in sizes,'Unknown chromosome in donor input'
                    dense=np.zeros(sizes[c],dtype=np.uint32)
                    # Source data concatenate eye/library blocks. A chromosome
                    # can recur; merge its earlier block within this run only.
                    if c in seen:
                        previous=np.load(dest/f'{c}.npz');dense[previous['positions']]=previous['counts']
                    seen.add(c);buffer=[]
                    for line in group:
                        _,a,b=line.split();a,b=int(a),int(b);left,right=a+delta_a,b+delta_b
                        assert 0<=left<sizes[c] and 0<=right<sizes[c]
                        buffer.extend((left,right));n+=1;endpoint_count+=2;merged.write(line)
                        peakinput.write(f'{c}\t{left}\t{left+1}\n{c}\t{right}\t{right+1}\n')
                        if len(buffer)>=200000:
                            np.add.at(dense,np.asarray(buffer,dtype=np.int64),1);buffer=[]
                    if buffer:np.add.at(dense,np.asarray(buffer,dtype=np.int64),1)
                    positions=np.flatnonzero(dense).astype(np.int32);values=dense[positions]
                    np.savez_compressed(dest/f'{c}.npz',positions=positions,counts=values);del dense
                    print(d,c,n,flush=True)
            for c in set(sizes)-seen:np.savez_compressed(dest/f'{c}.npz',positions=np.array([],dtype=np.int32),counts=np.array([],dtype=np.uint32))
            assert n==int(row['retained_fragment_records'])
            bwtotal=write_bigwig(out/f'{d}.insertions.bw',sizes,[dest]);assert bwtotal==2*n
            records[d]=dict(fragments=n,insertions=bwtotal,bigwig=str(out/f'{d}.insertions.bw'))
    Path(str(bed)+'.part').replace(bed);Path(str(fragments)+'.part').replace(fragments)
    pooled=write_bigwig(out/'muller.insertions.bw',sizes,sources);assert pooled==sum(r['insertions'] for r in records.values())
    report=dict(status='PASS',donors=records,pooled_insertions=pooled,pooled_fragments=pooled//2,
                pooled_bigwig=str(out/'muller.insertions.bw'),peak_input=str(bed),pipeline_fragment_input=str(fragments),
                fragment_input_coordinate_convention='Unmodified Cell Ranger fragment endpoints; upstream pipeline must detect shift. BigWigs are already corrected.',
                plus_endpoint_delta=delta_a,minus_endpoint_delta=delta_b,weighting='One molecule per row, two insertions per molecule; no read-support weighting or cross-cell deduplication.')
    Path('reports/preparation/tracks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':run()
