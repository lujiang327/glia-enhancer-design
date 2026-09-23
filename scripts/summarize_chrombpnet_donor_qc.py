#!/usr/bin/env python3
"""Aggregate donor-resolved ChromBPNet and matched-bias held-out metrics."""

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


DONORS = ("LGS1", "LGS2", "LGS3", "LVG1")
GROUPS = ("peaks", "nonpeaks", "peaks_and_nonpeaks")


def load_metrics(path):
    data = json.loads(path.read_text())
    result = {}
    for group in GROUPS:
        counts = data["counts_metrics"][group]
        profile = data["profile_metrics"][group]
        result[group] = {
            "counts_pearsonr": float(counts["pearsonr"]),
            "counts_spearmanr": float(counts["spearmanr"]),
            "counts_mse_uncalibrated": float(counts["mse"]),
            "median_jsd": float(profile["median_jsd"]),
            "median_normalized_jsd": float(profile["median_norm_jsd"]),
        }
    values = [value for group in result.values() for value in group.values()]
    if not all(math.isfinite(value) for value in values):
        raise ValueError("Non-finite donor metric in {}".format(path))
    return result


def load_donor_summary(path):
    with path.open(newline="") as handle:
        return {row["donor"]: row for row in csv.DictReader(handle, delimiter="\t")}


def describe(values):
    mean = statistics.mean(values)
    sd = statistics.stdev(values) if len(values) > 1 else 0.0
    return {
        "mean": mean,
        "standard_deviation": sd,
        "minimum": min(values),
        "maximum": max(values),
        "range": max(values) - min(values),
        "coefficient_of_variation": sd / abs(mean) if mean else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qc-root", type=Path, required=True)
    parser.add_argument("--donor-summary", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args()

    donor_info = load_donor_summary(args.donor_summary)
    records = []
    nested = {}
    cautions = []
    for donor in DONORS:
        donor_dir = args.qc_root / donor / "evaluation"
        model = load_metrics(donor_dir / "muller_chrombpnet_metrics.json")
        bias = load_metrics(donor_dir / "muller_scaled_bias_metrics.json")
        info = donor_info[donor]
        nested[donor] = {
            "cells": int(info["unambiguous_MG_cells"]),
            "retained_fragment_records": int(info["retained_fragment_records"]),
            "median_retained_fragments_per_cell": float(info["median_retained_fragments_per_cell"]),
            "pooled_TSS_proxy": float(info["pooled_TSS_proxy"]),
            "model_metrics": model,
            "matched_scaled_bias_metrics": bias,
            "improvement_over_scaled_bias": {},
        }
        for group in GROUPS:
            gains = {
                "counts_pearsonr_gain": model[group]["counts_pearsonr"] - bias[group]["counts_pearsonr"],
                "counts_spearmanr_gain": model[group]["counts_spearmanr"] - bias[group]["counts_spearmanr"],
                "median_jsd_reduction": bias[group]["median_jsd"] - model[group]["median_jsd"],
                "median_normalized_jsd_gain": model[group]["median_normalized_jsd"] - bias[group]["median_normalized_jsd"],
            }
            nested[donor]["improvement_over_scaled_bias"][group] = gains
            record = {
                "donor": donor,
                "group": group,
                "cells": nested[donor]["cells"],
                "retained_fragment_records": nested[donor]["retained_fragment_records"],
            }
            for prefix, values in (("model", model[group]), ("bias", bias[group])):
                for key, value in values.items():
                    record["{}_{}".format(prefix, key)] = value
            record.update(gains)
            records.append(record)

        peak_gain = nested[donor]["improvement_over_scaled_bias"]["peaks"]
        if model["peaks"]["counts_pearsonr"] <= 0:
            cautions.append("{} has non-positive peak count Pearson correlation.".format(donor))
        for key, value in peak_gain.items():
            if value <= 0:
                cautions.append("{} full model did not improve peak {} over bias.".format(donor, key))

    metric_specs = {
        "peak_counts_pearsonr": ("peaks", "counts_pearsonr"),
        "peak_counts_spearmanr": ("peaks", "counts_spearmanr"),
        "peak_median_jsd": ("peaks", "median_jsd"),
        "peak_median_normalized_jsd": ("peaks", "median_normalized_jsd"),
        "combined_counts_pearsonr": ("peaks_and_nonpeaks", "counts_pearsonr"),
        "combined_counts_spearmanr": ("peaks_and_nonpeaks", "counts_spearmanr"),
        "combined_median_jsd": ("peaks_and_nonpeaks", "median_jsd"),
        "combined_median_normalized_jsd": ("peaks_and_nonpeaks", "median_normalized_jsd"),
    }
    across = {}
    for name, (group, metric) in metric_specs.items():
        across[name] = describe([nested[d]["model_metrics"][group][metric] for d in DONORS])

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "status": "DONOR_RESULTS_READY_FOR_REVIEW",
        "fold": 0,
        "seed": 42,
        "negative_sampling_ratio": 1.0,
        "test_chromosomes": ["chr1", "chr3", "chr6"],
        "donors": nested,
        "across_donor_summary": across,
        "cautions": cautions,
        "count_mse_policy": (
            "Absolute count MSE is retained as an uncalibrated diagnostic but excluded "
            "from donor pass/fail because the checkpoint count scale reflects pooled depth."
        ),
        "scope": (
            "Donor-resolved held-out evaluation of a pooled model. This tests whether "
            "sequence-derived count ranking and profile shape generalize across donors; "
            "it is not leave-one-donor-out retraining."
        ),
    }
    args.output_json.write_text(json.dumps(summary, indent=2) + "\n")
    with args.output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
