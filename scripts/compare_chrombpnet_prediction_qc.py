#!/usr/bin/env python3
"""Compare held-out ChromBPNet metrics with a matched scaled-bias baseline."""
import argparse
import json
import math
from pathlib import Path


GROUPS = ("peaks", "nonpeaks", "peaks_and_nonpeaks")


def extract(metrics):
    counts = metrics["counts_metrics"]
    profiles = metrics["profile_metrics"]
    result = {}
    for group in GROUPS:
        if group not in counts or group not in profiles:
            raise ValueError(f"Missing prediction group: {group}")
        result[group] = {
            "counts_pearsonr": float(counts[group]["pearsonr"]),
            "counts_spearmanr": float(counts[group]["spearmanr"]),
            "counts_mse": float(counts[group]["mse"]),
            "median_jsd": float(profiles[group]["median_jsd"]),
            "median_normalized_jsd": float(profiles[group]["median_norm_jsd"]),
        }
    values = [value for group in result.values() for value in group.values()]
    if not all(math.isfinite(value) for value in values):
        raise ValueError("Non-finite held-out prediction metric")
    return result


def run(model_path, bias_path, output_path):
    model = extract(json.loads(Path(model_path).read_text()))
    bias = extract(json.loads(Path(bias_path).read_text()))

    # Positive deltas always mean that the full model improved over matched bias.
    deltas = {}
    for group in GROUPS:
        deltas[group] = {
            "counts_pearsonr_gain": (
                model[group]["counts_pearsonr"] - bias[group]["counts_pearsonr"]
            ),
            "counts_spearmanr_gain": (
                model[group]["counts_spearmanr"] - bias[group]["counts_spearmanr"]
            ),
            "counts_mse_reduction": (
                bias[group]["counts_mse"] - model[group]["counts_mse"]
            ),
            "median_jsd_reduction": (
                bias[group]["median_jsd"] - model[group]["median_jsd"]
            ),
            # Despite its upstream name, median_norm_jsd is a clipped
            # goodness score: 1 - JSD/JSD(observed, uniform). Higher is better.
            "median_normalized_jsd_gain": (
                model[group]["median_normalized_jsd"]
                - bias[group]["median_normalized_jsd"]
            ),
        }

    cautions = []
    peak = deltas["peaks"]
    for metric in (
        "counts_pearsonr_gain",
        "counts_spearmanr_gain",
        "counts_mse_reduction",
        "median_jsd_reduction",
        "median_normalized_jsd_gain",
    ):
        if peak[metric] <= 0:
            cautions.append(
                f"Full model did not improve peak {metric} over matched scaled bias."
            )

    # Background profile shape is expected to remain largely bias-driven, but
    # worse background count calibration can indicate that the biological
    # branch adds accessibility where it should be close to zero.
    nonpeak = deltas["nonpeaks"]
    for metric in (
        "counts_pearsonr_gain",
        "counts_spearmanr_gain",
        "counts_mse_reduction",
    ):
        if nonpeak[metric] <= 0:
            cautions.append(
                f"Full model did not improve nonpeak {metric} over matched scaled bias."
            )

    summary = {
        "status": "HELD_OUT_NUMERICAL_RESULTS_READY_FOR_REVIEW",
        "fold": 0,
        "test_chromosomes": ["chr1", "chr3", "chr6"],
        "model_metrics": model,
        "matched_scaled_bias_metrics": bias,
        "improvement_over_scaled_bias": deltas,
        "cautions": cautions,
        "interpretation": (
            "Positive improvement values favor the full model. This comparison "
            "uses identical held-out regions for both models and deliberately "
            "does not impose a universal correlation cutoff."
        ),
        "scope": (
            "Initial held-out aggregate numerical comparison only. Per-chromosome, "
            "depth-stratified, donor, seed, attribution, and motif QC remain required."
        ),
    }
    Path(output_path).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-metrics", required=True)
    parser.add_argument("--bias-metrics", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    run(args.model_metrics, args.bias_metrics, args.output)
