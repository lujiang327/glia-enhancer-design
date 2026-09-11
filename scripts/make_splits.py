#!/usr/bin/env python3
"""Freeze published retina chromosome partitions, restricted to chr1–22,X.

Source: doi:10.5281/zenodo.6796067, jobscripts/train/submit_training_fold*.sh.
These partitions are reused, not the legacy model or its performance claims.
"""
import json
from pathlib import Path


def run():
    pairs=[
        (['chr1','chr3','chr6'],['chr8','chr20']),
        (['chr2','chr8','chr9','chr16'],['chr12','chr17']),
        (['chr4','chr11','chr12','chr15'],['chr22','chr7']),
        (['chr5','chr10','chr14','chr18','chr20','chr22'],['chr6','chr21']),
        (['chr7','chr13','chr17','chr19','chr21','chrX'],['chr18','chr10']),
    ]
    out=Path('config/splits');out.mkdir(exist_ok=True)
    universe=[f'chr{i}' for i in range(1,23)]+['chrX']
    for i,(test,valid) in enumerate(pairs):
        train=[c for c in universe if c not in test+valid]
        (out/f'fold_{i}.json').write_text(json.dumps(dict(train=train,valid=valid,test=test),indent=2)+'\n')
    (out/'README.md').write_text('Five predeclared chromosome folds from the archived retina-models training scripts (https://doi.org/10.5281/zenodo.6796067). chrY was excluded to match the QC chromosome policy. Regenerate with `python3 scripts/make_splits.py`. The same split applies to peaks, negatives and bias fitting. Record interval hashes and inspect homologous regions separately. No model has been trained on these folds in this project.\n')


if __name__=='__main__':run()
