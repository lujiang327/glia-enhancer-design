#!/usr/bin/env python3
"""Build matched observed browser tracks from checksum-pinned donor pseudobulks."""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np
import pyBigWig
from prepare_browser_tracks import digest, bin_interval_counts


def count_endpoints(path,sizes,bin_bp,expected_fragments):
    counts={chrom:np.zeros((length+bin_bp-1)//bin_bp,dtype=np.uint64) for chrom,length in sizes.items()}
    current=None;buffer=[];records=0
    def flush():
        if buffer:
            np.add.at(counts[current],np.asarray(buffer,dtype=np.int64),1)
            buffer.clear()
    with gzip.open(path,'rt') as stream:
        for line in stream:
            fields=line.split()
            if len(fields)!=3: raise ValueError('Expected three-column donor BEDPE')
            chrom,start,end=fields[0],int(fields[1]),int(fields[2])
            # Match the verified MG track: endpoints at start and end, no extra shift.
            if chrom not in sizes or not 0<=start<end<sizes[chrom]:
                raise ValueError('Invalid endpoint under MG convention: '+line.strip())
            if current!=chrom:
                flush();current=chrom
            buffer.extend((start//bin_bp,end//bin_bp));records+=1
            if len(buffer)>=200000: flush()
            if records%5000000==0: print(path.name,records,'fragments',flush=True)
    flush()
    if records!=expected_fragments: raise ValueError('Fragment count mismatch')
    if sum(int(c.sum()) for c in counts.values())!=2*records: raise ValueError('Endpoint count mismatch')
    return counts


def write_track(path,counts,sizes,bin_bp,insertions):
    with pyBigWig.open(str(path)+'.part','w') as bw:
        bw.addHeader(list(sizes.items()))
        for chrom,length in sizes.items():
            for begin in range(0,len(counts[chrom]),50000):
                positions=np.arange(begin,min(begin+50000,len(counts[chrom])),dtype=np.int64)*bin_bp
                ends=np.minimum(positions+bin_bp,length)
                values=counts[chrom][begin:begin+len(positions)].astype(np.float64)/(ends-positions)*(1e6/insertions)
                bw.addEntries([chrom]*len(positions),positions,ends=ends,values=values)
    Path(str(path)+'.part').replace(path)
    with pyBigWig.open(str(path)) as check:
        if not np.isclose(check.header()['sumData'],1e6,rtol=2e-5): raise ValueError('Normalization integral failure')
    return dict(file=path.name,sha256=digest(path),insertions=insertions,bytes=path.stat().st_size)


def compare_mg_reference(counts,source,sizes,bin_bp):
    """Require exact integer-bin agreement with the previously validated MG track."""
    with pyBigWig.open(str(source)) as bw:
        if bw.chroms()!=sizes: raise ValueError('MG reference assembly mismatch')
        for chrom,length in sizes.items():
            for begin in range(0,length,bin_bp*50000):
                end=min(begin+bin_bp*50000,length)
                reference=bin_interval_counts(bw.intervals(chrom,begin,end),begin,end,bin_bp)
                observed=counts[chrom][begin//bin_bp:begin//bin_bp+len(reference)]
                if not np.array_equal(reference,observed): raise ValueError('MG reference endpoint/bin mismatch')
    return 'PASS_EXACT_ALL_GENOME_INTEGER_BIN_MATCH'


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--group',required=True)
    p.add_argument('--pseudobulk-dir',type=Path,required=True)
    p.add_argument('--chrom-sizes',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--mg-reference',type=Path)
    a=p.parse_args()
    config=json.loads(a.config.read_text());settings=config['groups'][a.group]
    if digest(a.chrom_sizes)!=config['chrom_sizes_sha256']: raise ValueError('Genome sizes hash mismatch')
    if a.group=='MG' and not a.mg_reference: p.error('MG requires --mg-reference')
    sizes=dict((c,int(n)) for c,n in (line.split() for line in a.chrom_sizes.read_text().splitlines()))
    bin_bp=config['bin_bp'];pooled={chrom:np.zeros((length+bin_bp-1)//bin_bp,dtype=np.uint64) for chrom,length in sizes.items()}
    a.output_dir.mkdir(parents=True,exist_ok=False)
    records=[];track_records=[];total_fragments=0
    for donor,expected in settings['donors'].items():
        path=a.pseudobulk_dir/(donor+'.fragments.bedpe.gz')
        if digest(path)!=expected['sha256']: raise ValueError('Pseudobulk checksum mismatch: '+donor)
        counts=count_endpoints(path,sizes,bin_bp,expected['fragments'])
        for chrom in sizes: pooled[chrom]+=counts[chrom]
        total_fragments+=expected['fragments']
        records.append(dict(donor=donor,**expected))
        if settings['keep_donor_tracks']:
            track_records.append(write_track(a.output_dir/(a.group+'_'+donor+'.IPM_per_base.100bp.bw'),counts,sizes,bin_bp,2*expected['fragments']))
        del counts
    reference_check=None
    if a.group=='MG':
        if digest(a.mg_reference)!=config['mg_reference_sha256']: raise ValueError('Raw MG reference hash mismatch')
        reference_check=compare_mg_reference(pooled,a.mg_reference,sizes,bin_bp)
    track_records.insert(0,write_track(a.output_dir/(a.group+'_pooled.IPM_per_base.100bp.bw'),pooled,sizes,bin_bp,2*total_fragments))
    summary=dict(status='OBSERVED_CELLTYPE_BROWSER_TRACKS_READY',group=a.group,assembly='hg38',bin_bp=bin_bp,
        fragments=total_fragments,insertions=2*total_fragments,donors=records,tracks=track_records,mg_reference_check=reference_check,
        config_sha256=digest(a.config),coordinate_policy='same source endpoints as validated MG tracks; no additional shift',
        unit='insertions per million retained insertions per base; bin mean; explicit zeros',
        limitations=['observed natural-locus ATAC, not reporter activity or mutant predictions',
            'library-size scaling does not estimate absolute accessibility per cell',
            'low-depth donors yield noisier visual peaks; retain original depth/cell information'])
    (a.output_dir/'browser_track_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__': main()
