#!/usr/bin/env python3
"""Empirical donor-specific Tn5 shift check using pinned ChromBPNet math.

Sample 100,000 fragment rows per donor, split into two 50,000-row checks.
Both ends of each sampled molecule are used. Sampling is seeded and uniform
within donor; rows with an N in either 40-bp context are excluded from the PWM.
"""
import csv
import json
from pathlib import Path
import numpy as np
from pyfaidx import Fasta
from vendor.shift_math import get_ref_pwms, compute_shift_ATAC


def pwm(sequences):
    letters=np.frombuffer(''.join(sequences).encode(),dtype='S1').reshape(-1,40)
    return np.stack([(letters==b).mean(axis=0) for b in [b'A',b'C',b'G',b'T']],axis=1)


def run():
    out=Path('reports/preparation');out.mkdir(exist_ok=True)
    refs=get_ref_pwms('config/ATAC.ref.motifs.txt');results={}
    with Fasta('data/raw/reference/hg38.canonical.fa',as_raw=True) as genome:
        for row in csv.DictReader(open('reports/qc/donor_summary.tsv'),delimiter='\t'):
            donor=row['donor'];total=int(row['retained_fragment_records']);rng=np.random.default_rng(42)
            wanted=set(rng.choice(total,size=min(100000,total),replace=False).tolist());chosen=[]
            with open(f'data/intermediate/donor_pseudobulks/{donor}.bedpe') as f:
                for i,line in enumerate(f):
                    if i in wanted:
                        c,a,b=line.split();chosen.append((c,int(a),int(b)))
                assert i+1==total
            rng.shuffle(chosen);checks=[]
            for subset in [chosen[:len(chosen)//2],chosen[len(chosen)//2:]]:
                plus=[];minus=[]
                for c,a,b in subset:
                    if a<20 or b+20>len(genome[c]):continue
                    left=genome[c][a-20:a+20].upper();right=genome[c][b-20:b+20].upper()
                    if set(left+right)<=set('ACGT'):plus.append(left);minus.append(right)
                assert len(plus)>10000
                pp,mp=pwm(plus),pwm(minus);shift=compute_shift_ATAC(*refs,pp,mp)
                checks.append(dict(n_fragments=len(plus),plus_shift=int(shift[0]),minus_shift=int(shift[1])))
                np.savez_compressed(out/f'{donor}.shift_check_{len(checks)}.npz',plus_pwm=pp,minus_pwm=mp)
            assert len({(r['plus_shift'],r['minus_shift']) for r in checks})==1,'Unstable shift estimates'
            results[donor]=checks;print(donor,checks,flush=True)
    unique={(r['plus_shift'],r['minus_shift']) for checks in results.values() for r in checks}
    assert len(unique)==1,'Donors have different coordinate shifts: process separately'
    plus,minus=unique.pop()
    record=dict(status='PASS',checks=results,seed=42,current_plus_shift=plus,current_minus_shift=minus,target_plus_shift=4,target_minus_shift=-4,
        plus_endpoint_delta=4-plus,minus_endpoint_delta=-4-minus-1,
        convention='Insertion coordinates: start+(4-observed_plus), end+(-4-observed_minus)-1; matches upstream bedtools genomecov -5.',policy=__doc__)
    (out/'fragment_shift.json').write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':run()
