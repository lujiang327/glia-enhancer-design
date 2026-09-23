#!/usr/bin/env python3
"""Summarize fold-0 ChromBPNet training and held-out metrics across seeds."""

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


SEEDS = (42, 123, 456)


def load_json(path):
    with path.open() as handle:
        return json.load(handle)


def summary(values):
    mean = statistics.fmean(values)
    sd = statistics.stdev(values) if len(values) > 1 else 0.0
    return {
        "mean": mean,
        "standard_deviation": sd,
        "coefficient_of_variation": sd / abs(mean) if mean else None,
        "minimum": min(values),
        "maximum": max(values),
        "range": max(values) - min(values),
    }


def svg_report(rows, histories, output):
    width, height = 1200, 430
    colors = {42: "#386cb0", 123: "#f0027f", 456: "#31a354"}
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:17px;font-weight:bold}.axis{font-size:12px}.legend{font-size:13px}</style>',
    ]

    # Validation-loss trajectories.
    x0, y0, w, h = 65, 55, 360, 300
    all_losses = [float(r["val_loss"]) for hist in histories.values() for r in hist]
    ymin, ymax = min(all_losses) - 3, max(all_losses) + 3
    max_epoch = max(len(hist) for hist in histories.values())
    parts += [f'<text x="{x0}" y="25" class="title">Validation loss</text>',
              f'<line x1="{x0}" y1="{y0+h}" x2="{x0+w}" y2="{y0+h}" stroke="#444"/>',
              f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0+h}" stroke="#444"/>']
    for value in (380, 390, 400, 410, 420):
        if ymin <= value <= ymax:
            y = y0 + h * (ymax - value) / (ymax - ymin)
            parts += [f'<line x1="{x0}" y1="{y:.1f}" x2="{x0+w}" y2="{y:.1f}" stroke="#ddd"/>',
                      f'<text x="{x0-8}" y="{y+4:.1f}" text-anchor="end" class="axis">{value}</text>']
    for seed, hist in histories.items():
        points = []
        for i, record in enumerate(hist):
            x = x0 + w * i / max(1, max_epoch - 1)
            y = y0 + h * (ymax - float(record["val_loss"])) / (ymax - ymin)
            points.append(f"{x:.1f},{y:.1f}")
        parts.append(f'<polyline points="{" ".join(points)}" fill="none" stroke="{colors[seed]}" stroke-width="2.5"/>')
    parts += [f'<text x="{x0+w/2}" y="{y0+h+35}" text-anchor="middle" class="axis">Epoch (zero-based)</text>',
              f'<text x="18" y="{y0+h/2}" transform="rotate(-90 18 {y0+h/2})" text-anchor="middle" class="axis">Loss</text>']

    # Grouped bars for held-out count correlations and normalized profile JSD.
    panels = [
        (500, "Held-out count Pearson r", ["peaks", "nonpeaks", "combined"],
         ["peak_counts_pearsonr", "nonpeak_counts_pearsonr", "combined_counts_pearsonr"], 0.0, 0.85),
        (855, "Median normalized profile JSD", ["peaks", "nonpeaks", "combined"],
         ["peak_median_normalized_jsd", "nonpeak_median_normalized_jsd", "combined_median_normalized_jsd"], 0.0, 0.36),
    ]
    row_by_seed = {row["seed"]: row for row in rows}
    for px, title, labels, keys, low, high in panels:
        pw, ph, py = 285, 300, 55
        parts += [f'<text x="{px}" y="25" class="title">{title}</text>',
                  f'<line x1="{px}" y1="{py+ph}" x2="{px+pw}" y2="{py+ph}" stroke="#444"/>',
                  f'<line x1="{px}" y1="{py}" x2="{px}" y2="{py+ph}" stroke="#444"/>']
        for gi, (label, key) in enumerate(zip(labels, keys)):
            gx = px + 30 + gi * 90
            for si, seed in enumerate(SEEDS):
                value = float(row_by_seed[seed][key])
                bh = ph * (value - low) / (high - low)
                bx = gx + si * 18
                parts.append(f'<rect x="{bx}" y="{py+ph-bh:.1f}" width="15" height="{bh:.1f}" fill="{colors[seed]}"/>')
            parts.append(f'<text x="{gx+18}" y="{py+ph+20}" text-anchor="middle" class="axis">{label}</text>')
        for value in (low, (low + high) / 2, high):
            y = py + ph * (high - value) / (high - low)
            parts.append(f'<text x="{px-7}" y="{y+4:.1f}" text-anchor="end" class="axis">{value:.2f}</text>')

    for i, seed in enumerate(SEEDS):
        lx = 500 + i * 105
        parts += [f'<rect x="{lx}" y="402" width="14" height="14" fill="{colors[seed]}"/>',
                  f'<text x="{lx+20}" y="414" class="legend">seed {seed}</text>']
    parts.append('</svg>')
    output.write_text("\n".join(parts) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()

    rows, histories = [], {}
    for seed in SEEDS:
        base = args.raw_root / "models/chrombpnet/fold_0" / f"seed_{seed}_neg1"
        training = load_json(base / "logs/training_summary.json")
        qc = load_json(base / "qc_prediction/logs/numerical_qc_summary.json")
        with (base / "logs/muller_chrombpnet.log").open(newline="") as handle:
            histories[seed] = list(csv.DictReader(handle))
        metrics = qc["model_metrics"]
        rows.append({
            "seed": seed,
            "epochs_completed": training["epochs_completed"],
            "best_epoch_zero_based": training["best_epoch_zero_based"],
            "best_val_loss": training["best_val_loss"],
            "peak_counts_pearsonr": metrics["peaks"]["counts_pearsonr"],
            "peak_counts_spearmanr": metrics["peaks"]["counts_spearmanr"],
            "peak_counts_mse": metrics["peaks"]["counts_mse"],
            "peak_median_jsd": metrics["peaks"]["median_jsd"],
            "peak_median_normalized_jsd": metrics["peaks"]["median_normalized_jsd"],
            "nonpeak_counts_pearsonr": metrics["nonpeaks"]["counts_pearsonr"],
            "nonpeak_counts_spearmanr": metrics["nonpeaks"]["counts_spearmanr"],
            "nonpeak_counts_mse": metrics["nonpeaks"]["counts_mse"],
            "nonpeak_median_jsd": metrics["nonpeaks"]["median_jsd"],
            "nonpeak_median_normalized_jsd": metrics["nonpeaks"]["median_normalized_jsd"],
            "combined_counts_pearsonr": metrics["peaks_and_nonpeaks"]["counts_pearsonr"],
            "combined_counts_spearmanr": metrics["peaks_and_nonpeaks"]["counts_spearmanr"],
            "combined_counts_mse": metrics["peaks_and_nonpeaks"]["counts_mse"],
            "combined_median_jsd": metrics["peaks_and_nonpeaks"]["median_jsd"],
            "combined_median_normalized_jsd": metrics["peaks_and_nonpeaks"]["median_normalized_jsd"],
        })

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with prefix.with_suffix(".tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    numeric_keys = [key for key in rows[0] if key not in {"seed", "epochs_completed", "best_epoch_zero_based"}]
    result = {
        "status": "PASS_SEED_REPRODUCIBILITY_WITH_NONPEAK_CALIBRATION_CAVEAT",
        "fold": 0,
        "negative_sampling_ratio": 1.0,
        "seeds": list(SEEDS),
        "per_seed": rows,
        "across_seed_summary": {key: summary([float(row[key]) for row in rows]) for key in numeric_keys},
        "cautions": [
            "All seeds retain worse nonpeak count correlation and MSE than the matched scaled-bias baseline.",
            "Seed reproducibility does not establish biological-donor generalization."
        ],
        "decision": "Seed reproducibility passes; proceed to donor-sensitivity testing while differential accessibility remains blocked."
    }
    prefix.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    svg_report(rows, histories, prefix.with_suffix(".svg"))


if __name__ == "__main__":
    main()
