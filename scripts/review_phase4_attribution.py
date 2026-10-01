#!/usr/bin/env python3
"""Apply the Phase 4 orientation-agreement gate and make compact QC outputs."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


AGREEMENT_MIN = 0.5
SIGN_CONCORDANCE_MIN = 0.8


def failure_reasons(row):
    reasons = []
    for head in ("count", "profile"):
        if row[f"{head}_attribution_parent_pearson"] < AGREEMENT_MIN:
            reasons.append(f"{head}_pearson_lt_{AGREEMENT_MIN}")
        if row[f"{head}_attribution_parent_cosine"] < AGREEMENT_MIN:
            reasons.append(f"{head}_cosine_lt_{AGREEMENT_MIN}")
        if row[f"{head}_attribution_parent_top10pct_sign_concordance"] < SIGN_CONCORDANCE_MIN:
            reasons.append(f"{head}_top10pct_sign_lt_{SIGN_CONCORDANCE_MIN}")
    if row["model_input_invalid_bases"] != 0:
        reasons.append("ambiguous_model_input_bases")
    return ";".join(reasons)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--agreement", required=True, type=Path)
    parser.add_argument("--output-table", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    parser.add_argument("--output-figure", required=True, type=Path)
    args = parser.parse_args()

    table = pd.read_csv(args.agreement, sep="\t")
    if len(table) != 1000 or table["peak_id"].nunique() != 1000:
        raise ValueError("Expected 1,000 unique candidate rows")
    table["attribution_qc_failure_reasons"] = table.apply(failure_reasons, axis=1)
    table["attribution_qc_pass"] = table["attribution_qc_failure_reasons"].eq("")

    args.output_table.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_table, sep="\t", index=False, compression="gzip")

    failure_counts = {}
    for value in table.loc[~table["attribution_qc_pass"], "attribution_qc_failure_reasons"]:
        for reason in value.split(";"):
            failure_counts[reason] = failure_counts.get(reason, 0) + 1
    class_counts = (
        table.groupby("genomic_class")["attribution_qc_pass"]
        .agg(total="size", passed="sum")
        .reset_index()
    )
    class_counts["failed"] = class_counts["total"] - class_counts["passed"]
    summary = {
        "status": "PASS_FOR_DUAL_ORIENTATION_PERTURBATION_WITH_FILTERING",
        "candidates": int(len(table)),
        "attribution_qc_pass": int(table["attribution_qc_pass"].sum()),
        "attribution_qc_fail": int((~table["attribution_qc_pass"]).sum()),
        "gate": {
            "count_parent_pearson_min": AGREEMENT_MIN,
            "count_parent_cosine_min": AGREEMENT_MIN,
            "count_parent_top10pct_sign_concordance_min": SIGN_CONCORDANCE_MIN,
            "profile_parent_pearson_min": AGREEMENT_MIN,
            "profile_parent_cosine_min": AGREEMENT_MIN,
            "profile_parent_top10pct_sign_concordance_min": SIGN_CONCORDANCE_MIN,
            "model_input_invalid_bases": 0,
        },
        "failure_counts_nonexclusive": dict(sorted(failure_counts.items())),
        "pass_by_genomic_class": class_counts.to_dict(orient="records"),
        "parent_count_orientation_flagged_total": int(table["parent_count_orientation_flag"].sum()),
        "parent_count_orientation_flagged_among_pass": int(
            (table["parent_count_orientation_flag"] & table["attribution_qc_pass"]).sum()
        ),
        "policy": (
            "Use only passing candidates for sequence perturbation; score every parent and variant "
            "in both orientations and require concordant improvement."
        ),
        "next_gate": "RUN_ATTRIBUTION_GUIDED_PERTURBATION_PILOT",
    }
    args.output_summary.write_text(json.dumps(summary, indent=2) + "\n")

    passed = table["attribution_qc_pass"].to_numpy()
    colors = np.where(passed, "#4c78a8", "#e45756")
    figure, axes = plt.subplots(2, 2, figsize=(11, 9))
    axes[0, 0].scatter(
        table["count_attribution_parent_pearson"],
        table["profile_attribution_parent_pearson"],
        c=colors, s=13, alpha=0.7,
    )
    axes[0, 0].axvline(AGREEMENT_MIN, color="black", linestyle="--", linewidth=1)
    axes[0, 0].axhline(AGREEMENT_MIN, color="black", linestyle="--", linewidth=1)
    axes[0, 0].set_xlabel("Count attribution F/RC Pearson")
    axes[0, 0].set_ylabel("Profile attribution F/RC Pearson")
    axes[0, 0].set_title("Parent attribution agreement")

    bins = np.linspace(-0.7, 1, 45)
    axes[0, 1].hist(table["count_attribution_parent_pearson"], bins=bins, alpha=0.65, label="count")
    axes[0, 1].hist(table["profile_attribution_parent_pearson"], bins=bins, alpha=0.65, label="profile")
    axes[0, 1].axvline(AGREEMENT_MIN, color="black", linestyle="--", linewidth=1)
    axes[0, 1].set_xlabel("Forward/aligned-RC Pearson")
    axes[0, 1].set_ylabel("Candidates")
    axes[0, 1].set_title("Agreement distributions")
    axes[0, 1].legend(frameon=False)

    axes[1, 0].scatter(
        table["candidate_rank"], table["profile_attribution_parent_pearson"],
        c=colors, s=13, alpha=0.7,
    )
    axes[1, 0].axhline(AGREEMENT_MIN, color="black", linestyle="--", linewidth=1)
    axes[1, 0].set_xlabel("Candidate rank")
    axes[1, 0].set_ylabel("Profile attribution F/RC Pearson")
    axes[1, 0].set_title("Profile agreement across ranking")

    class_counts = class_counts.sort_values("passed", ascending=False)
    axes[1, 1].bar(class_counts["genomic_class"], class_counts["passed"], label="pass", color="#4c78a8")
    axes[1, 1].bar(
        class_counts["genomic_class"], class_counts["failed"],
        bottom=class_counts["passed"], label="fail", color="#e45756",
    )
    axes[1, 1].set_ylabel("Candidates")
    axes[1, 1].set_title("QC outcome by genomic class")
    axes[1, 1].tick_params(axis="x", rotation=20)
    axes[1, 1].legend(frameon=False)

    figure.tight_layout()
    figure.savefig(args.output_figure, dpi=180)
    plt.close(figure)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
