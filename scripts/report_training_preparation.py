#!/usr/bin/env python3
"""Plot measured GC matching and summarize preparation without opening model gate."""
import csv
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', str(Path('data/intermediate/matplotlib').resolve()))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def run():
    root = Path('reports/preparation')
    regions = json.loads((root/'regions.json').read_text())
    tracks = json.loads((root/'tracks.json').read_text())
    validation = json.loads((root/'input_validation.json').read_text())
    gc_review = json.loads((root/'gc_review.json').read_text())
    assert gc_review['parameters'] == regions['parameters'], 'Renew GC review after parameter changes'
    fig, axes = plt.subplots(5, 3, figsize=(12, 13), sharex=True, sharey=True)
    table = ['| Fold / split | Peaks | Negatives | Exact GC bucket | Mean GC gap | Maximum gap |',
             '|---|---:|---:|---:|---:|---:|']
    for i, (fold, record) in enumerate(regions['folds'].items()):
        with open(f'data/intermediate/training/{fold}.gc_matches.tsv') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        for j, split in enumerate(['train', 'valid', 'test']):
            values = [r for r in rows if r['split'] == split]
            ax = axes[i, j]
            for field, label in [('positive_gc_bucket', 'Peak'), ('negative_gc_bucket', 'Background')]:
                hist = np.bincount([int(r[field]) for r in values], minlength=101)
                ax.plot(np.arange(101)/100, hist/hist.sum(), label=label, alpha=.85)
            ax.set_title(f'{fold} / {split}'); ax.set_xlim(.15, .85)
            if i == 4: ax.set_xlabel('GC fraction (rounded to 0.01)')
            if j == 0: ax.set_ylabel('Fraction of windows')
            m = record['metrics'][split]
            table.append(f"| {fold} / {split} | {m['peaks']:,} | {m['negatives']:,} | {m['exact_bucket_fraction']:.1%} | {m['mean_gc_gap']:.4f} | {m['max_gc_gap']:.2f} |")
    axes[0,0].legend(); fig.tight_layout()
    for ext in ['png', 'svg']: fig.savefig(root/f'gc_matching.{ext}', dpi=160)
    plt.close(fig)
    counts = regions['peak_filter_counts']
    message = f'''# Müller glia training-input preparation

Input-integrity validation: **{validation['status']}**. Model status: **NOT TRAINED**.
The GC review accepts these backgrounds for an exploratory pilot, with residual
mismatches recorded below. Proceed only if input-integrity validation also passes.
These results do not establish a reliable Müller glia model. Differential
accessibility and sequence design remain blocked.

- Four donor tracks; {tracks['pooled_fragments']:,} retained fragment records and
  {tracks['pooled_insertions']:,} pooled insertions.
- hg38 analysis-set reference checksum verified; chr1–22 and X retained.
- Both independent shift checks for every donor recover native +4/−5 endpoints.
- {counts['raw_peaks']:,} relaxed genome-wide peaks before filtering;
  {counts['retained_peaks']:,} training peaks after sequence/context filtering.
- {regions['candidate_background_count']:,} eligible background grid windows.
- Five folds with disjoint train/validation/test chromosomes; regions checked for ACGT,
  chromosome bounds, blacklist overlap and full 500-bp jitter context.
- Exact central input sequence/reverse-complement leakage is checked across
  splits. Near-homology and offset homology are not exhaustively evaluated.

## GC matching

GC gaps are absolute fractions (0.01 = one percentage point). Backgrounds are
selected without replacement within each split, using the nearest available
bucket when exact matching is exhausted. The figure shows every fold and split;
peak distributions are repeated to match the requested background ratio.

The initial 1,000-bp grid was rejected because of high-GC mismatch; its outputs
are archived. The final 100-bp grid gives at least
{gc_review['minimum_exact_bucket_fraction']:.1%} exact bucket matches per split.
At most {gc_review['maximum_fraction_gap_gt_0_05']:.2%} of pairs per split have a
gap greater than 0.05. Individual gaps still reach
{gc_review['maximum_individual_gc_gap']:.2f}; GC-stratified held-out performance
must be inspected. Acceptance is a project judgment for a pilot, not a published
ChromBPNet threshold. Neighboring backgrounds can overlap within a split.

''' + '\n'.join(table) + '''

![GC distributions](gc_matching.png)

## Reproducibility and remaining work

See [the preparation guide](../../docs/training_preparation.md) for commands,
coordinate conventions, MACS3 and GC-sampling adaptations, and transfer steps.
Parameters, shift evidence, counts, validation and SHA256 manifests accompany
this report. Raw fragments and derived large files are retained outside Git.

QC exclusions and donor caveats remain those in the
[Phase 1 report](../phase1_qc_report.md). Input quality does not resolve potential
RNA contamination or donor imbalance. Training still requires a verified Linux
GPU runtime, bias fitting, multiple seeds, held-out profile/count evaluation,
donor generalization checks and attribution/motif review. No checkpoint exists.
'''
    (root/'README.md').write_text(message)


if __name__ == '__main__':
    run()
