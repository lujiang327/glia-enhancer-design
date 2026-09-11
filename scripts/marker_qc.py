#!/usr/bin/env python3
"""Check RNA markers only on unambiguously mapped published MG barcodes.

Does not cluster, transfer labels, relabel cells, or perform differential tests.
"""
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path


def run(root=Path('.')):
    raw=root/'data/raw';out=root/'reports/qc'
    samples={r['rna_prefix']:r for r in csv.DictReader(open(root/'config/samples.tsv'),delimiter='\t')}
    audit={r['bare_barcode']:r for r in csv.DictReader(open(out/'barcode_audit.tsv'),delimiter='\t')}
    with gzip.open(raw/'GSM5866081_barcodes.tsv.gz','rt') as f: cells=[l.strip() for l in f]
    markers=['GLUL','RLBP1','SLC1A3','SOX9','AQP4','GFAP','PDE6A','RHO','ARR3','GRM6','GAD1','C1QA']
    with gzip.open(raw/'GSM5866081_features.tsv.gz','rt') as f:
        features=[l.strip().split('\t') for l in f]
    indices={i+1:r[1] for i,r in enumerate(features) if r[1] in markers}
    totals=[0]*len(cells);values={m:{} for m in markers}
    with gzip.open(raw/'GSM5866081_matrix.mtx.gz','rt') as f:
        line=next(f)
        if not line.startswith('%%MatrixMarket matrix coordinate'):raise ValueError('Expected coordinate MTX')
        line=next(f)
        while line.startswith('%'):line=next(f)
        nfeat,ncells,nnz=map(int,line.split())
        assert (nfeat,ncells)==(len(features),len(cells))
        seen=0
        for line in f:
            i,j,v=map(int,line.split());j-=1;totals[j]+=v;seen+=1
            if i in indices:values[indices[i]][j]=v
        assert seen==nnz
    groups=defaultdict(list);linked=[]
    for j,cell in enumerate(cells):
        sample,barcode=cell.split('_',1);r=audit.get(barcode)
        if r and r['n_matching_RNA_samples']=='1':
            assert sample==r['matching_RNA_sample_prefixes']
            donor=samples[sample]['donor'];groups['MG_'+donor].append(j);groups['MG_all'].append(j)
            linked.append([cell,barcode,sample,donor,samples[sample]['eye'],'Muller glia','published_fragment_file_unique_RNA_match'])
        elif r:groups['ambiguous_excluded'].append(j)
        else:groups['non_MG_barcode_background'].append(j)
    with open(out/'cells.unambiguous.tsv','w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['cell_id','bare_barcode','rna_prefix','donor','eye','published_label','label_evidence']);w.writerows(linked)
    with open(out/'marker_summary.tsv','w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['group','gene','n_cells','fraction_detected','mean_log1p_CP10k'])
        for group,js in sorted(groups.items()):
            for m in markers:
                v=values[m]
                w.writerow([group,m,len(js),sum(v.get(j,0)>0 for j in js)/len(js),sum(math.log1p(v.get(j,0)*10000/totals[j]) if totals[j] else 0 for j in js)/len(js)])
    (out/'marker_qc.json').write_text(json.dumps(dict(rna_cells=len(cells),features=len(features),nnz=nnz,groups={g:len(js) for g,js in groups.items()},missing_markers=sorted(set(markers)-set(indices.values())),policy=__doc__),indent=2)+'\n')
    print((out/'marker_qc.json').read_text())


if __name__=='__main__':run()
