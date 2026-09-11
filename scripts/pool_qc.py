#!/usr/bin/env python3
"""Pooled hg38 QC; writes BEDPE for peak calling, no differential testing.

Assumes input is a deduplicated 10x fragment file with already shifted ends.
Retains chr1-22,X; excludes full fragments overlapping blacklist intervals.
TSS proxy: strand-oriented start and end-1 insertions, unique RefGene TSSs,
mean signal +/-50 bp divided by mean outer 100 bp within +/-2 kb.
This is not the ArchR per-cell TSS enrichment statistic.
"""
import argparse
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path


def merged_bed(path):
    intervals = defaultdict(list)
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith(('#', 'track')):
                continue
            c, a, b, *_ = line.split()
            intervals[c].append((int(a), int(b)))
    result = {}
    for c, values in intervals.items():
        merged = []
        for a, b in sorted(values):
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        result[c] = ([a for a,b in merged], [b for a,b in merged])
    return result


def overlaps(intervals, chrom, start, end):
    starts, ends = intervals.get(chrom, ([], []))
    i = bisect_left(starts, end) - 1
    return i >= 0 and ends[i] > start


def run(args):
    args.out.mkdir(parents=True, exist_ok=True)
    canonical = {f'chr{i}' for i in range(1,23)} | {'chrX'}
    sizes = dict((c,int(n)) for c,n in (l.split() for l in args.sizes.read_text().splitlines()))
    blacklist = merged_bed(args.blacklist)
    tss = defaultdict(set)
    with gzip.open(args.refgene, 'rt') as f:
        for line in f:
            row = line.split('\t')
            c, strand = row[2:4]
            pos = int(row[4]) if strand == '+' else int(row[5]) - 1
            if c in canonical and 2000 <= pos < sizes[c]-2000 and not overlaps(blacklist,c,pos-2000,pos+2001):
                tss[(c,strand)].add(pos)
    tss = {k:sorted(v) for k,v in tss.items()}
    profile = [0]*4001
    counts = Counter()
    per_barcode = Counter()
    with gzip.open(args.fragments,'rt') as f, (args.out/'muller.filtered.bedpe').open('w') as bed:
        for line in f:
            if line.startswith('#'): continue
            c,a,b,barcode,_ = line.split()
            a,b=int(a),int(b)
            counts['input'] += 1
            if c not in canonical:
                counts['noncanonical'] += 1
                continue
            if a < 0 or b > sizes[c] or b <= a:
                counts['invalid_coordinates'] += 1
                continue
            if overlaps(blacklist,c,a,b):
                counts['blacklist_overlap'] += 1
                continue
            counts['retained'] += 1
            per_barcode[barcode] += 1
            bed.write(f'{c}\t{a}\t{b}\n')
            for strand, sign in [('+',1),('-',-1)]:
                positions = tss.get((c,strand), [])
                for endpoint in (a,b-1):
                    lo=bisect_left(positions,endpoint-2000)
                    hi=bisect_right(positions,endpoint+2000)
                    for p in positions[lo:hi]:
                        profile[(endpoint-p)*sign+2000] += 1
            if counts['input'] % 10000000 == 0: print(dict(counts),flush=True)
    flank = (sum(profile[:100])+sum(profile[-100:]))/200
    center = sum(profile[1950:2051])/101
    result = dict(counts=counts, tss_count=sum(map(len,tss.values())),
                  pooled_tss_enrichment_proxy=center/flank if flank else None,
                  tss_definition=__doc__, coordinate_convention='10x shifted start and end-1; no further shift')
    (args.out/'pool_qc.json').write_text(json.dumps(result,indent=2)+'\n')
    for filename, header, rows in [
        ('tss_profile.tsv',['offset','insertions'],enumerate(profile,-2000)),
        ('filtered_barcode_counts.tsv',['bare_barcode','retained_records'],sorted(per_barcode.items()))]:
        with (args.out/filename).open('w') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(header);w.writerows(rows)
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['fragments','blacklist','refgene','sizes','out']:
        p.add_argument('--'+key,type=Path,required=True)
    run(p.parse_args())
