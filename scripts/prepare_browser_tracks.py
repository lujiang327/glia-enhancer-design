#!/usr/bin/env python3
"""Export compact observed MG browser tracks without changing training inputs."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import pyBigWig


def digest(path):
    result=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):
            result.update(block)
    return result.hexdigest()


def bin_interval_counts(intervals, start, end, bin_bp):
    """Integrate clipped intervals into fixed bins, including boundary crossings."""
    count = (end-start+bin_bp-1)//bin_bp
    result = np.zeros(count, dtype=np.float64)
    if not intervals:
        return result
    values = np.asarray(intervals, dtype=np.float64)
    left = np.maximum(values[:,0].astype(np.int64),start)-start
    right = np.minimum(values[:,1].astype(np.int64),end)-start
    v = values[:,2]
    valid = right>left
    left,right,v = left[valid],right[valid],v[valid]
    first,last = left//bin_bp,(right-1)//bin_bp
    np.add.at(result,first,v*(np.minimum(right,(first+1)*bin_bp)-left))
    cross = last>first
    np.add.at(result,last[cross],v[cross]*(right[cross]-last[cross]*bin_bp))
    diff = np.zeros(count+1,dtype=np.float64)
    np.add.at(diff,first[cross]+1,v[cross]*bin_bp)
    np.add.at(diff,last[cross],-v[cross]*bin_bp)
    result += np.cumsum(diff[:-1])
    return result


def export(source, destination, expected_insertions, bin_bp):
    total=0.
    with pyBigWig.open(str(source)) as src, pyBigWig.open(str(destination)+'.part','w') as dest:
        sizes=src.chroms()
        if not np.isclose(src.header()['sumData'],expected_insertions,rtol=1e-10):
            raise ValueError('Raw track insertion sum mismatch: '+str(source))
        dest.addHeader(list(sizes.items()))
        for chrom,length in sizes.items():
            for start in range(0,length,bin_bp*50000):
                end=min(start+bin_bp*50000,length)
                intervals=src.intervals(chrom,start,end)
                counts=bin_interval_counts(intervals,start,end,bin_bp)
                total += counts.sum()
                positions=np.arange(start,end,bin_bp,dtype=np.int64)
                ends=np.minimum(positions+bin_bp,length)
                normalized=counts/(ends-positions)*(1e6/expected_insertions)
                # Include zero bins explicitly so the viewer sees zero, not unknown.
                dest.addEntries([chrom]*len(positions),positions,ends=ends,values=normalized)
            print(destination.name,chrom,flush=True)
    if not np.isclose(total,expected_insertions,rtol=1e-10):
        raise ValueError('Binned integration did not conserve insertion counts')
    Path(str(destination)+'.part').replace(destination)
    with pyBigWig.open(str(destination)) as check:
        if not np.isclose(check.header()['sumData'],1e6,rtol=2e-5):
            raise ValueError('Normalized signal integral mismatch')
        return dict(source=str(source),file=destination.name,source_insertions=int(total),
            source_sha256=digest(source),
            normalized_signal_integral=check.header()['sumData'],
            sha256=digest(destination),bytes=destination.stat().st_size)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--bin-bp',type=int,default=100)
    a=p.parse_args()
    if a.bin_bp<=0: p.error('Bin width must be positive')
    a.output_dir.mkdir(parents=True,exist_ok=False)
    source=Path('reports/preparation/tracks.json')
    config=json.loads(source.read_text())
    tracks=[('MG_pooled',Path(config['pooled_bigwig']),config['pooled_insertions'])]
    tracks.extend(('MG_'+name,Path(data['bigwig']),data['insertions']) for name,data in config['donors'].items())
    records=[]
    for name,path,total in tracks:
        records.append(export(path,a.output_dir/(name+'.IPM_per_base.'+str(a.bin_bp)+'bp.bw'),total,a.bin_bp))
    parents=pd.read_csv('reports/phase4/parent_scoring/parent_predictions.tsv.gz',sep='\t').sort_values('candidate_rank')
    if len(parents)!=1000 or parents.peak_id.duplicated().any(): raise ValueError('Parent BED integrity failure')
    bed=parents[['chrom','start','end']].copy()
    bed['name']=['rank'+str(rank).zfill(4)+'__'+pid for rank,pid in zip(parents.candidate_rank,parents.peak_id)]
    bed['score']=0
    bed['strand']='.'
    bed.sort_values(['chrom','start']).to_csv(a.output_dir/'MG_top1000_parents.hg38.bed',sep='\t',index=False,header=False)
    parents[['candidate_rank','peak_id','chrom','start','end','genomic_class','nearest_gene']].to_csv(a.output_dir/'candidate_index.tsv',sep='\t',index=False)
    first=parents.iloc[0]
    locus=f'{first.chrom}:{max(1,int(first.start)-1000+1)}-{int(first.end)+1000}'
    (a.output_dir/'first_locus.txt').write_text(locus+'\n')
    summary=dict(status='OBSERVED_MG_BROWSER_TRACKS_READY_OFFTARGET_TRACKS_NOT_INCLUDED',assembly='hg38',
        signal='observed ATAC insertion density, not model prediction',unit='insertions per million retained insertions per base',
        bin_bp=a.bin_bp,normalization='bin insertion count / actual bin width * 1e6 / whole-track retained insertion count',
        display_note='100bp bin mean for overview; raw single-base tracks remain available at their source paths',
        shifts='inherited from previously validated training tracks; no additional endpoint shift',
        zero_bins='explicitly included',tracks=records,candidate_regions=len(parents),first_locus=locus,
        limitations=['MG tracks alone do not demonstrate off-target specificity','donor tracks use separate depth normalization',
            '100bp bins do not resolve single-base footprints','designed variants have no observed ATAC signal',
            'global library normalization does not equal absolute accessibility per cell'],
        input_summary_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    (a.output_dir/'browser_track_manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
    shutil.copy2('data/raw/reference/hg38.canonical.chrom.sizes',a.output_dir/'hg38.canonical.chrom.sizes')
    for file in sorted(a.output_dir.iterdir()):
        print(file.name,file.stat().st_size,flush=True)
    (a.output_dir/'checksums.sha256').write_text(''.join(digest(file)+'  '+file.name+'\n' for file in sorted(a.output_dir.iterdir()) if file.name!='checksums.sha256'))

if __name__=='__main__': main()
