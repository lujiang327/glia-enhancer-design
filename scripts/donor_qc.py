#!/usr/bin/env python3
"""Build identifiable-donor pseudobulks from the unambiguous barcode subset.

50-kb fragment-midpoint bins are QC features, not accessible peaks or DA tests.
"""
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
import statistics
from pool_qc import merged_bed, overlaps


def run(root=Path('.')):
    out=root/'reports/qc';raw=root/'data/raw';dest=root/'data/intermediate/donor_pseudobulks';dest.mkdir(parents=True,exist_ok=True)
    samples={r['rna_prefix']:r for r in csv.DictReader(open(root/'config/samples.tsv'),delimiter='\t')}
    audit=list(csv.DictReader(open(out/'barcode_audit.tsv'),delimiter='\t'))
    mapping={r['bare_barcode']:samples[r['matching_RNA_sample_prefixes']]['donor'] for r in audit if r['n_matching_RNA_samples']=='1'}
    donor_names=sorted(set(mapping.values()))
    canonical={f'chr{i}' for i in range(1,23)}|{'chrX'}
    sizes=dict((c,int(n)) for c,n in (l.split() for l in (raw/'hg38.chrom.sizes').read_text().splitlines()))
    blacklist=merged_bed(raw/'hg38-blacklist.v2.bed.gz')
    tss=defaultdict(set)
    with gzip.open(raw/'hg38.refGene.txt.gz','rt') as f:
        for line in f:
            r=line.split('\t');c,strand=r[2:4];p=int(r[4]) if strand=='+' else int(r[5])-1
            if c in canonical and 2000<=p<sizes[c]-2000 and not overlaps(blacklist,c,p-2000,p+2001):tss[(c,strand)].add(p)
    tss={k:sorted(v) for k,v in tss.items()}
    profiles={d:[0]*4001 for d in donor_names};bins={d:Counter() for d in donor_names};cells=Counter();flow=Counter()
    handles={d:open(dest/f'{d}.bedpe','w') for d in donor_names}
    try:
        with gzip.open(raw/'GSM5866073_Mullerglia_frags.tsv.gz','rt') as f:
            for line in f:
                if line.startswith('#'):continue
                c,a,b,barcode,_=line.split();a,b=int(a),int(b);flow['input']+=1
                if barcode not in mapping:flow['ambiguous_excluded']+=1;continue
                if c not in canonical or a<0 or b>sizes[c] or b<=a or overlaps(blacklist,c,a,b):flow['genomic_filter_excluded']+=1;continue
                d=mapping[barcode];flow['retained']+=1;cells[barcode]+=1;bins[d][(c,((a+b)//2)//50000)]+=1
                handles[d].write(f'{c}\t{a}\t{b}\n')
                for strand,sign in [('+',1),('-',-1)]:
                    positions=tss.get((c,strand),[])
                    for endpoint in (a,b-1):
                        lo=bisect_left(positions,endpoint-2000);hi=bisect_right(positions,endpoint+2000)
                        for p in positions[lo:hi]:profiles[d][(endpoint-p)*sign+2000]+=1
                if flow['input']%10000000==0:print(dict(flow),flush=True)
    finally:
        for handle in handles.values():handle.close()
    with open(out/'donor_summary.tsv','w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['donor','unambiguous_MG_cells','retained_fragment_records','median_retained_fragments_per_cell','pooled_TSS_proxy','scope'])
        for d in donor_names:
            vals=[cells[b] for b,donor in mapping.items() if donor==d];p=profiles[d];flank=(sum(p[:100])+sum(p[-100:]))/200
            w.writerow([d,len(vals),sum(vals),statistics.median(vals),(sum(p[1950:2051])/101)/flank if flank else '', 'unambiguous_subset_only'])
    with open(out/'donor_50kb_bins.tsv','w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['chrom','start',*donor_names])
        for c in sorted(canonical):
            for b in range((sizes[c]+49999)//50000):
                if not overlaps(blacklist,c,b*50000,min((b+1)*50000,sizes[c])):
                    w.writerow([c,b*50000,*[bins[d][(c,b)] for d in donor_names]])
    with open(out/'donor_tss_profiles.tsv','w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['offset',*donor_names]);w.writerows([i-2000,*[profiles[d][i] for d in donor_names]] for i in range(4001))
    (out/'donor_qc.json').write_text(json.dumps(dict(flow=flow,bin_size=50000,tss_count=sum(map(len,tss.values())),policy=__doc__),indent=2)+'\n')
    print((out/'donor_summary.tsv').read_text())


if __name__=='__main__':run()
