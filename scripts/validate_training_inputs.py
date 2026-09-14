#!/usr/bin/env python3
"""Independently validate prepared sequences, tracks and chromosome partitions.

This gate concerns input integrity only. It cannot validate a trained model.
"""
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import pyBigWig
from pyfaidx import Fasta
from pool_qc import merged_bed,overlaps
from prepare_reference import digest


def read_regions(path):
    with open(path) as f:
        for line in f:
            r=line.split();assert len(r)==10
            yield r[0],int(r[1])+int(r[9]),r


def run():
    out=Path('reports/preparation');config=json.loads(Path('config/preparation.json').read_text());half=config['input_length']//2;margin=half+config['max_jitter']
    ref=json.loads((out/'reference.json').read_text());tracks=json.loads((out/'tracks.json').read_text());regions=json.loads((out/'regions.json').read_text())
    shift=json.loads((out/'fragment_shift.json').read_text());assert shift['status']=='PASS'
    fasta=Path(ref['fasta']);assert digest(fasta)==ref['fasta_sha256']
    sizes=dict((c,int(n)) for c,n in (l.split() for l in Path('data/raw/reference/hg38.canonical.chrom.sizes').read_text().splitlines()))
    blacklist=merged_bed('data/raw/hg38-blacklist.v2.bed.gz');rawpeaks=merged_bed('data/intermediate/training/peaks/muller_relaxed_peaks.narrowPeak')
    donors={d:pyBigWig.open(r['bigwig']) for d,r in tracks['donors'].items()};pooled=pyBigWig.open(tracks['pooled_bigwig'])
    assert pooled.chroms()==sizes
    for d,bw in donors.items():
        assert bw.chroms()==sizes
        assert int(round(bw.header()['sumData']))==2*tracks['donors'][d]['fragments']
    assert int(round(pooled.header()['sumData']))==tracks['pooled_insertions']==2*tracks['pooled_fragments']
    rng=np.random.default_rng(7123)
    for c,size in sizes.items():
        for a in rng.integers(0,size-1000,size=10):
            a=int(a);expected=sum(np.nan_to_num(bw.values(c,a,a+1000,numpy=True)) for bw in donors.values())
            observed=np.nan_to_num(pooled.values(c,a,a+1000,numpy=True));assert np.array_equal(observed,expected)
    for bw in donors.values():bw.close()
    pooled.close()
    peaklist=list(read_regions(regions['peaks']));assert len(peaklist)==regions['peak_filter_counts']['retained_peaks']
    reverse=str.maketrans('ACGT','TGCA');fold_reports={}
    with Fasta(str(fasta),as_raw=True) as genome:
        peak_hashes=[]
        for c,center,r in peaklist:
            assert center-margin>=0 and center+margin<=sizes[c]
            assert not overlaps(blacklist,c,center-margin,center+margin)
            context=genome[c][center-margin:center+margin].upper();assert len(context)==2*margin and set(context)<=set('ACGT')
            seq=context[config['max_jitter']:config['max_jitter']+config['input_length']]
            peak_hashes.append(hashlib.sha256(min(seq,seq.translate(reverse)[::-1]).encode()).digest())
        for fold,record in regions['folds'].items():
            split=json.loads(Path(f'config/splits/{fold}.json').read_text());mapping={c:s for s,cs in split.items() for c in cs}
            assert sum(map(len,split.values()))==len(mapping)==len(sizes) and set(mapping)==set(sizes)
            seen_sequences={};cross_split_duplicates=[]
            for (c,center,r),key in zip(peaklist,peak_hashes):
                sp=mapping[c]
                if key in seen_sequences and seen_sequences[key]!=sp:cross_split_duplicates.append((c,center,'peak'))
                seen_sequences[key]=sp
            totals=Counter();seen=set()
            for c,center,r in read_regions(record['negatives']):
                sp=mapping[c];totals[sp]+=1;assert (c,center) not in seen;seen.add((c,center))
                assert center-margin>=0 and center+margin<=sizes[c]
                assert not overlaps(blacklist,c,center-margin,center+margin)
                assert not overlaps(rawpeaks,c,center-margin,center+margin)
                context=genome[c][center-margin:center+margin].upper();assert len(context)==2*margin and set(context)<=set('ACGT')
                seq=context[config['max_jitter']:config['max_jitter']+config['input_length']]
                key=hashlib.sha256(min(seq,seq.translate(reverse)[::-1]).encode()).digest()
                if key in seen_sequences and seen_sequences[key]!=sp:cross_split_duplicates.append((c,center,'negative'))
                seen_sequences[key]=sp
            for sp in split:
                n=sum(mapping[c]==sp for c,center,r in peaklist)
                ratio=config['background_ratio_test'] if sp=='test' else config['background_ratio_train_valid']
                assert totals[sp]==n*ratio==record['metrics'][sp]['negatives']
            fold_reports[fold]=dict(region_counts=dict(totals),exact_sequence_or_reverse_complement_cross_split_duplicates=len(cross_split_duplicates),duplicate_examples=cross_split_duplicates[:10])
            print(fold,fold_reports[fold],flush=True)
    duplicate_total=sum(r['exact_sequence_or_reverse_complement_cross_split_duplicates'] for r in fold_reports.values())
    validation=dict(status='PASS' if duplicate_total==0 else 'REVIEW_SEQUENCE_DUPLICATES',
        scope='Input integrity, not model reliability',reference_sha256_verified=True,track_totals_verified=True,
        pooled_equals_donor_sum_at_230_sampled_windows=True,all_regions_valid_acgt_with_jitter=True,
        chromosome_partitions_disjoint=True,folds=fold_reports,
        limitations=['Near-homology and jitter-offset homologous sequences are not exhaustively assessed','GC matching deviations require report review','Bias fitting and all model validation remain pending'])
    (out/'input_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    assert duplicate_total==0,'Exact sequence leakage found; resolve before training'
    paths=[str(fasta),str(fasta)+'.fai','data/raw/reference/hg38.canonical.chrom.sizes','data/raw/hg38-blacklist.v2.bed.gz',tracks['pipeline_fragment_input'],tracks['pooled_bigwig'],regions['peaks']]
    paths += [r['bigwig'] for r in tracks['donors'].values()]
    paths += [r['negatives'] for r in regions['folds'].values()]
    # Include donor fragments for leave-one-donor-out fits without recovering cell IDs.
    paths += [str(p) for p in sorted(Path('data/intermediate/donor_pseudobulks').glob('*.bedpe'))]
    manifest=[dict(path=p,bytes=Path(p).stat().st_size,sha256=digest(p)) for p in paths]
    (out/'training_input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    Path('config/training_transfer_files.txt').write_text('\n'.join(paths)+'\n')
    plan=json.loads(Path('config/training_plan.json').read_text())
    plan.update(reference_fasta=str(fasta),chrom_sizes='data/raw/reference/hg38.canonical.chrom.sizes',training_peak_set=regions['peaks'],gc_matched_nonpeak_set={f:r['negatives'] for f,r in regions['folds'].items()},
        insertion_bigwig=tracks['pooled_bigwig'],input_fragment_file=tracks['pipeline_fragment_input'],max_jitter=config['max_jitter'],status='INPUTS_PREPARED_MODEL_NOT_TRAINED')
    Path('config/training_plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    print('Input validation complete; model gate remains closed.')


if __name__=='__main__':run()
