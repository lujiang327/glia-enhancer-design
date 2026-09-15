#!/usr/bin/env python3
"""Apply the documented ChromBPNet numerical screen to bias prediction metrics."""
import argparse
import json
import math
from pathlib import Path


def run(metrics_path, output_path):
    metrics = json.loads(Path(metrics_path).read_text())
    counts = metrics['counts_metrics']
    profiles = metrics['profile_metrics']
    required_groups = ('nonpeaks', 'peaks')
    for group in required_groups:
        if group not in counts or group not in profiles:
            raise ValueError(f'Missing prediction group: {group}')

    values = {
        'nonpeak_counts_pearsonr': float(counts['nonpeaks']['pearsonr']),
        'nonpeak_counts_spearmanr': float(counts['nonpeaks']['spearmanr']),
        'nonpeak_counts_mse': float(counts['nonpeaks']['mse']),
        'nonpeak_median_jsd': float(profiles['nonpeaks']['median_jsd']),
        'nonpeak_median_normalized_jsd': float(profiles['nonpeaks']['median_norm_jsd']),
        'peak_counts_pearsonr': float(counts['peaks']['pearsonr']),
        'peak_counts_spearmanr': float(counts['peaks']['spearmanr']),
        'peak_counts_mse': float(counts['peaks']['mse']),
        'peak_median_jsd': float(profiles['peaks']['median_jsd']),
        'peak_median_normalized_jsd': float(profiles['peaks']['median_norm_jsd']),
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError('Non-finite bias prediction metric')

    failures = []
    cautions = []
    if values['nonpeak_counts_pearsonr'] <= 0:
        failures.append('Held-out nonpeak count Pearson correlation is not positive.')
    if values['peak_counts_pearsonr'] <= -0.5:
        failures.append('Peak count Pearson correlation is at or below the upstream rejection boundary (-0.5).')
    elif values['peak_counts_pearsonr'] <= -0.3:
        cautions.append('Peak count Pearson correlation is in the upstream caution range (-0.5, -0.3].')

    if failures:
        status = 'FAIL_NUMERICAL_BIAS_QC'
    elif cautions:
        status = 'CAUTION_NUMERICAL_QC_MOTIF_QC_REQUIRED'
    else:
        status = 'PASS_NUMERICAL_QC_MOTIF_QC_PENDING'

    summary = {
        'status': status,
        'fold': 0,
        'test_chromosomes': ['chr1', 'chr3', 'chr6'],
        'metrics': values,
        'failures': failures,
        'cautions': cautions,
        'criteria': {
            'nonpeak_counts_pearsonr': '> 0',
            'peak_counts_pearsonr_preferred': '> -0.3',
            'peak_counts_pearsonr_reject': '<= -0.5',
        },
        'scope': 'Held-out numerical screen only. Attribution and motif-leakage QC remain required before bias-model acceptance.'
    }
    Path(output_path).write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metrics', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    run(args.metrics, args.output)
