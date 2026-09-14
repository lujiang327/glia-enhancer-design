#!/usr/bin/env python3
"""Filter relaxed peaks and sample fold-specific GC-matched genomic negatives.

Native Mac implementation of the ChromBPNet GC-bucket strategy, with bounded
nearest-bucket fallback instead of the upstream unbounded random walk. Adds
ACGT and full-jitter context checks. All matching discrepancies are retained.
"""
from bisect import bisect_left
from collections import Counter,defaultdict
import csv
import json
from pathlib import Path
import numpy as np
from pyfaidx import Fasta
from pool_qc import merged_bed,overlaps


def sequence_ok(seq,a,b):
    return 0<=a<b<=len(seq) and seq.count('N',a,b)==0


def gc_bucket(seq,a,b):
    return int(round(100*round((seq.count('G',a,b)+seq.count('C',a,b))/(b-a),2)))


def choose_bucket(available,wanted):
    if not available:raise ValueError('No unused GC-matching candidates remain in this chromosome split')
    return min(available,key=lambda x:(abs(x-wanted),x))


def run():
    config=json.loads(Path('config/preparation.json').read_text());half=config['input_length']//2;margin=half+config['max_jitter'];stride=config['background_stride']
    out=Path('data/intermediate/training');reportdir=Path('reports/preparation')
    peakpath=out/'peaks/muller_relaxed_peaks.narrowPeak'
    rawlines=defaultdict(list)
    for line in peakpath.read_text().splitlines():rawlines[line.split('\t')[0]].append(line.split('\t'))
    blacklist=merged_bed('data/raw/hg38-blacklist.v2.bed.gz');rawpeaks=merged_bed(peakpath)
    sizes=dict((c,int(n)) for c,n in (l.split() for l in Path('data/raw/reference/hg38.canonical.chrom.sizes').read_text().splitlines()))
    positives=[];candidates=[];flow=Counter();chrom_index={c:i for i,c in enumerate(sizes)}
    with Fasta('data/raw/reference/hg38.canonical.fa',as_raw=True) as genome:
        for c in sizes:
            seq=genome[c][:].upper();seen=set()
            for r in rawlines[c]:
                flow['raw_peaks']+=1;a,b=int(r[1]),int(r[2]);summit=a+int(r[9])
                if not a<=summit<b:flow['invalid_summit']+=1;continue
                if summit in seen:flow['duplicate_summit']+=1;continue
                if not sequence_ok(seq,summit-margin,summit+margin):flow['boundary_or_N_context']+=1;continue
                if overlaps(blacklist,c,max(0,a-margin),min(len(seq),b+margin)):flow['blacklist_context']+=1;continue
                seen.add(summit);bucket=gc_bucket(seq,summit-half,summit+half);positives.append((r,bucket));flow['retained_peaks']+=1
            before=len(candidates)
            for a in range(0,len(seq)-config['input_length']+1,stride):
                summit=a+half;lo,hi=summit-margin,summit+margin
                if not sequence_ok(seq,lo,hi):continue
                if overlaps(blacklist,c,lo,hi) or overlaps(rawpeaks,c,lo,hi):continue
                candidates.append((chrom_index[c],summit,gc_bucket(seq,a,a+2*half)))
            print(c,'peaks',len(positives),'background candidates',len(candidates)-before,flush=True)
    assert positives and candidates
    filtered=out/'muller.training_peaks.narrowPeak'
    filtered.write_text(''.join('\t'.join(r)+'\n' for r,_ in positives))
    with (out/'peak_gc.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['chrom','summit','gc_bucket']);w.writerows((r[0],int(r[1])+int(r[9]),g) for r,g in positives)
    candidate_array=np.array(candidates,dtype=np.int32);del candidates
    np.savez_compressed(out/'background_candidates.npz',candidates=candidate_array,chromosomes=np.array(list(sizes)))
    folds={};names=list(sizes)
    for foldpath in sorted(Path('config/splits').glob('fold_*.json')):
        split=json.loads(foldpath.read_text());rng=np.random.default_rng(config['background_seed']+int(foldpath.stem.split('_')[1]));mapping={c:s for s,cs in split.items() for c in cs}
        assert set(mapping)==set(sizes);bins=defaultdict(list)
        for i,(ci,center,gc) in enumerate(candidate_array):bins[(mapping[names[ci]],int(gc))].append(i)
        for key in bins:rng.shuffle(bins[key])
        available={s:{g for (sp,g),vs in bins.items() if sp==s and vs} for s in split}
        selected=[];metrics={};matches=[]
        for r,wanted in positives:
            sp=mapping[r[0]];ratio=config['background_ratio_test'] if sp=='test' else config['background_ratio_train_valid']
            for _ in range(ratio):
                bucket=choose_bucket(available[sp],wanted);index=bins[(sp,bucket)].pop()
                if not bins[(sp,bucket)]:available[sp].remove(bucket)
                ci,center,gc=map(int,candidate_array[index]);selected.append((ci,center,index));matches.append((sp,wanted,gc))
        assert len({i for _,_,i in selected})==len(selected)
        negative=out/f'{foldpath.stem}.negatives.bed'
        with negative.open('w') as f:
            for k,(ci,center,index) in enumerate(sorted(selected)):
                f.write(f'{names[ci]}\t{center-half}\t{center+half}\tbg_{index}\t0\t.\t0\t0\t0\t{half}\n')
        with (out/f'{foldpath.stem}.gc_matches.tsv').open('w') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['split','positive_gc_bucket','negative_gc_bucket']);w.writerows(matches)
        for sp in split:
            foreground=np.array([a for s,a,b in matches if s==sp]);background=np.array([b for s,a,b in matches if s==sp]);diff=np.abs(foreground-background)/100
            metrics[sp]=dict(peaks=sum(mapping[r[0]]==sp for r,g in positives),negatives=len(diff),exact_bucket_fraction=float(np.mean(diff==0)),mean_gc_gap=float(diff.mean()),max_gc_gap=float(diff.max()),fraction_gap_gt_0_05=float(np.mean(diff>.05)),gc_histogram_total_variation=float(np.abs(np.bincount(foreground,minlength=101)-np.bincount(background,minlength=101)).sum()/2/len(diff)))
        folds[foldpath.stem]=dict(negatives=str(negative),metrics=metrics);print(foldpath.stem,metrics,flush=True)
    report=dict(status='GENERATED_REQUIRES_VALIDATION',peak_filter_counts=dict(flow),candidate_background_count=len(candidate_array),peaks=str(filtered),folds=folds,parameters=config,method=__doc__)
    (reportdir/'regions.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':run()
