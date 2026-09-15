#!/usr/bin/env python3
"""Plot ChromBPNet fold-specific 1 kb peak and nonpeak count distributions."""
import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np
import pyBigWig

matplotlib.use('Agg')
import matplotlib.pyplot as plt


QUANTILES = (0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 1)


def read_regions(path, chromosomes):
    regions = []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip().split('\t')
            if fields[0] in chromosomes:
                center = int(fields[1]) + int(float(fields[9]))
                regions.append((fields[0], center - 500, center + 500))
    return regions


def window_counts(bigwig, regions):
    values = np.empty(len(regions), dtype=np.float64)
    for index, (chromosome, start, end) in enumerate(regions):
        value = bigwig.stats(chromosome, start, end, type='sum', exact=True)[0]
        values[index] = 0 if value is None else value
    return values


def describe(values):
    estimates = np.quantile(values, QUANTILES)
    return {
        'windows': int(values.size),
        'zero_fraction': float(np.mean(values == 0)),
        'mean': float(np.mean(values)),
        'quantiles': {str(q): float(v) for q, v in zip(QUANTILES, estimates)},
    }


def run(args):
    split = json.loads(Path(args.fold).read_text())
    chromosomes = set(split['train'] + split['valid'])
    peak_regions = read_regions(args.peaks, chromosomes)
    nonpeak_regions = read_regions(args.nonpeaks, chromosomes)

    bigwig = pyBigWig.open(args.bigwig)
    try:
        peak_counts = window_counts(bigwig, peak_regions)
        nonpeak_counts = window_counts(bigwig, nonpeak_regions)
    finally:
        bigwig.close()

    count_cutoff = float(np.quantile(peak_counts, 0.01) * args.bias_threshold)
    below_cutoff = nonpeak_counts[nonpeak_counts < count_cutoff]
    lower = float(np.quantile(below_cutoff, 1 - args.outlier_threshold))
    upper = float(np.quantile(below_cutoff, args.outlier_threshold))
    retained = below_cutoff[(below_cutoff > lower) & (below_cutoff < upper)]

    result = {
        'scope': 'Fold 0 train and validation chromosomes; test chromosomes excluded.',
        'window_width_bp': 1000,
        'peak_counts': describe(peak_counts),
        'nonpeak_counts': describe(nonpeak_counts),
        'bias_filter': {
            'factor_times_peak_1pct': args.bias_threshold,
            'upper_count_cutoff_strict': count_cutoff,
            'windows_below_cutoff': int(below_cutoff.size),
            'lower_outlier_bound_strict': lower,
            'upper_outlier_bound_strict': upper,
            'windows_retained': int(retained.size),
            'retained_median': float(np.median(retained)),
            'raw_counts_loss_weight': float(np.median(retained) / 10),
        },
    }
    Path(args.output_json).write_text(json.dumps(result, indent=2) + '\n')

    colors = {'Peaks': '#0072B2', 'Nonpeaks': '#D55E00'}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)

    maximum = max(np.quantile(peak_counts, 0.999), np.quantile(nonpeak_counts, 0.999))
    bins = np.linspace(0, np.log10(maximum + 1), 55)
    for label, values in [('Peaks', peak_counts), ('Nonpeaks', nonpeak_counts)]:
        axes[0].hist(
            np.log10(values + 1), bins=bins, histtype='step', linewidth=2,
            weights=np.full(values.size, 100 / values.size), color=colors[label], label=label
        )
    tick_counts = np.array([0, 1, 3, 10, 30, 100, 300, 1000, 3000])
    tick_counts = tick_counts[tick_counts <= maximum]
    axes[0].set_xticks(np.log10(tick_counts + 1), tick_counts)
    axes[0].set_xlabel('Tn5 insertions per 1 kb window')
    axes[0].set_ylabel('Windows per bin (%)')
    axes[0].set_title('Training + validation distributions')
    axes[0].legend(frameon=False)

    exact = np.arange(0, 51)
    percentages = np.array([np.mean(nonpeak_counts == value) * 100 for value in exact])
    axes[1].bar(exact, percentages, width=0.9, color=colors['Nonpeaks'])
    axes[1].axvline(count_cutoff, color='black', linestyle='--', linewidth=1.5,
                    label=f'upper cutoff = {count_cutoff:g}')
    axes[1].axvline(np.median(retained), color='#009E73', linestyle=':', linewidth=2,
                    label=f'retained median = {np.median(retained):g}')
    axes[1].set_xlim(-0.5, 50.5)
    axes[1].set_xlabel('Tn5 insertions per 1 kb nonpeak window')
    axes[1].set_ylabel('All nonpeak windows (%)')
    axes[1].set_title('Bias-region count range')
    axes[1].legend(frameon=False)

    fig.suptitle('Müller glia ChromBPNet bias-input counts (fold 0)', fontsize=13)
    for path in (args.output_png, args.output_svg):
        fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bigwig', default='data/intermediate/training/muller.insertions.bw')
    parser.add_argument('--peaks', default='data/intermediate/training/muller.training_peaks.narrowPeak')
    parser.add_argument('--nonpeaks', default='data/intermediate/training/fold_0.negatives.bed')
    parser.add_argument('--fold', default='config/splits/fold_0.json')
    parser.add_argument('--bias-threshold', type=float, default=0.5)
    parser.add_argument('--outlier-threshold', type=float, default=0.9999)
    parser.add_argument('--output-json', default='reports/preparation/bias_count_distribution_fold0.json')
    parser.add_argument('--output-png', default='reports/preparation/bias_count_distribution_fold0.png')
    parser.add_argument('--output-svg', default='reports/preparation/bias_count_distribution_fold0.svg')
    run(parser.parse_args())
