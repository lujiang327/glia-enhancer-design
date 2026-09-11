#!/usr/bin/env python3
"""FRiP on merged diagnostic MACS3 intervals, counting each fragment once."""
import csv
import json
from pathlib import Path
from pool_qc import merged_bed, overlaps


def run():
    out=Path('reports/qc');peaks=out/'peaks/muller_diagnostic_peaks.narrowPeak'
    blacklist=merged_bed('data/raw/hg38-blacklist.v2.bed.gz')
    kept=[]
    for line in peaks.read_text().splitlines():
        c,a,b,*_=line.split();a,b=int(a),int(b)
        if not overlaps(blacklist,c,a,b):kept.append(line)
    filtered=out/'peaks/muller_diagnostic.no_blacklist.narrowPeak';filtered.write_text('\n'.join(kept)+'\n')
    intervals=merged_bed(filtered)
    paths={'full_pool':out/'muller.filtered.bedpe'}
    paths.update({p.stem:p for p in Path('data/intermediate/donor_pseudobulks').glob('*.bedpe')})
    rows=[]
    for name,path in paths.items():
        total=inside=0
        with path.open() as f:
            for line in f:
                c,a,b=line.split();a,b=int(a),int(b);total+=1
                if overlaps(intervals,c,a,b):inside+=1
        rows.append([name,total,inside,inside/total if total else None]);print(name,total,inside,flush=True)
    with (out/'frip.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['scope','retained_fragment_records','fragments_overlapping_peak','FRiP']);w.writerows(rows)
    (out/'peak_qc.json').write_text(json.dumps(dict(diagnostic_peaks_before_blacklist=sum(1 for _ in peaks.open()),diagnostic_peaks_after_blacklist=len(kept),merged_peak_intervals=sum(len(v[0]) for v in intervals.values()),method='MACS3 BEDPE q=0.01; hs; keep-dup all. Full published MG pool. FRiP = any fragment overlap / retained records. Not final training peaks or an independent validation peak set.'),indent=2)+'\n')


if __name__=='__main__':run()
