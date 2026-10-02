#!/usr/bin/env python3
"""Verify full-scan edits, pilot reproducibility, and provisional edit coverage."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from run_phase4_perturbation_pilot import edit_screen_pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--pilot", required=True, type=Path)
    parser.add_argument("--parents", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    output = args.results_dir
    table = pd.read_csv(output / "perturbation_full.tsv.gz", sep="\t")
    coverage = pd.read_csv(output / "perturbation_full_parent_coverage.tsv.gz", sep="\t")
    pilot = pd.read_csv(args.pilot, sep="\t")
    parents = pd.read_csv(args.parents, sep="\t")
    config = json.loads(args.config.read_text())
    assert len(coverage) == config["parent_count"] == 947
    assert coverage.peak_id.nunique() == 947
    assert len(table) == table.variant_id.nunique() == 9470
    assert set(table.peak_id) == set(coverage.peak_id)
    assert np.isfinite(table.select_dtypes("number").to_numpy()).all()
    for row in table.itertuples():
        assert len(row.parent_sequence) == len(row.variant_sequence) == 500
        changed = [i for i, (a, b) in enumerate(zip(row.parent_sequence, row.variant_sequence)) if a != b]
        assert changed == [row.parent_position_0based]
        assert row.parent_sequence[changed[0]] == row.ref and row.variant_sequence[changed[0]] == row.alt
        assert row.genomic_position_0based == row.start + changed[0]
        for field in ("parent", "variant"):
            assert hashlib.sha256(getattr(row, field + "_sequence").encode()).hexdigest() == getattr(row, field + "_sequence_sha256")
    screen = config["edit_screen"]
    passing = edit_screen_pass(table.robust_effect_lower_bound,
                              table.forward_profile_jsd_vs_parent, table.reverse_profile_jsd_vs_parent,
                              screen["minimum_intended_logcount_effect_each_orientation"],
                              screen["maximum_profile_jsd_vs_parent_each_orientation"])
    assert np.array_equal(passing.to_numpy(), table.edit_screen_pass.to_numpy())
    for kind in ("gain", "loss_control"):
        count = table[(table.design_class == kind) & passing].groupby("peak_id").size()
        assert np.array_equal(coverage.peak_id.map(count).fillna(0).astype(int), coverage[kind + "_screen_pass"])

    repeated = pilot.merge(table, on="variant_id", suffixes=("_pilot", "_full"), validate="one_to_one")
    assert len(repeated) == len(pilot) == 1000
    max_differences = {}
    for field in ("delta_forward_logcount", "delta_reverse_complement_logcount", "delta_mean_logcount",
                  "forward_profile_jsd_vs_parent", "reverse_profile_jsd_vs_parent"):
        difference = float(abs(repeated[field + "_pilot"] - repeated[field + "_full"]).max())
        assert difference <= 1e-6
        max_differences[field] = difference

    selected = table[passing].sort_values(
        ["candidate_rank", "design_class", "robust_effect_lower_bound", "variant_id"],
        ascending=[True, True, False, True],
    ).drop_duplicates(["peak_id", "design_class"])
    selected.to_csv(output / "best_passing_single_edits.tsv.gz", sep="\t", index=False, compression="gzip")
    exceptions = coverage[(coverage.gain_screen_pass == 0) | (coverage.loss_control_screen_pass == 0)]
    exceptions = exceptions.merge(parents, on=["peak_id", "candidate_rank"], validate="one_to_one")
    exceptions.to_csv(output / "parents_without_passing_edit.tsv.gz", sep="\t", index=False, compression="gzip")

    summary = {
        "status": "PASS_FOR_PROVISIONAL_EDIT_PRIORITIZATION_LIBRARY_ANNOTATION_PENDING",
        "parents": len(coverage), "variants": len(table),
        "parents_with_passing_gain": int((coverage.gain_screen_pass > 0).sum()),
        "parents_with_passing_loss_control": int((coverage.loss_control_screen_pass > 0).sum()),
        "parents_with_both": int(((coverage.gain_screen_pass > 0) & (coverage.loss_control_screen_pass > 0)).sum()),
        "parents_without_gain": int((coverage.gain_screen_pass == 0).sum()),
        "parents_without_loss_control": int((coverage.loss_control_screen_pass == 0).sum()),
        "best_passing_edit_rows": len(selected),
        "repeated_pilot_variants": len(repeated), "pilot_max_absolute_differences": max_differences,
        "by_design_class": {},
        "next_gate": "ANNOTATE_REPEAT_OVERLAP_AND_MAPPABILITY_BEFORE_LIBRARY_SELECTION",
    }
    for kind, group in table.groupby("design_class"):
        max_jsd = screen["maximum_profile_jsd_vs_parent_each_orientation"]
        summary["by_design_class"][kind] = {
            "variants": len(group), "screen_pass": int(group.edit_screen_pass.sum()),
            "median_delta_mean_logcount": float(group.delta_mean_logcount.median()),
            "nomination_direct_effect_spearman": float(group.nomination_count_hypothetical_delta.corr(group.delta_mean_logcount, method="spearman")),
            "profile_jsd_failures": int(((group.forward_profile_jsd_vs_parent > max_jsd) | (group.reverse_profile_jsd_vs_parent > max_jsd)).sum()),
        }
    (output / "full_review_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for kind, color in (("gain", "#4c78a8"), ("loss_control", "#e45756")):
        group = table[table.design_class == kind]
        axes[0, 0].scatter(group.delta_forward_logcount, group.delta_reverse_complement_logcount, s=5, alpha=.25, color=color, label=kind)
        values = coverage[kind + "_screen_pass"].value_counts().reindex(range(6), fill_value=0)
        offset = -.18 if kind == "gain" else .18
        axes[0, 1].bar(np.arange(6) + offset, values, width=.36, color=color, label=kind)
        axes[1, 0].hist(group.robust_effect_lower_bound, bins=50, alpha=.6, color=color, label=kind)
        axes[1, 1].scatter(group.candidate_rank, group.robust_effect_lower_bound, s=5, alpha=.2, color=color)
    for axis in (axes[0, 0],):
        axis.axhline(0, linewidth=.8, color="black")
        axis.axvline(0, linewidth=.8, color="black")
    axes[0, 0].set(xlabel="Forward count effect (log units)", ylabel="RC count effect (log units)", title="Direct effect agreement")
    axes[0, 0].legend(frameon=False)
    axes[0, 1].set(xlabel="Passing edits per parent", ylabel="Parents", title="Passing-edit coverage", xticks=range(6))
    axes[0, 1].legend(frameon=False)
    min_effect = screen["minimum_intended_logcount_effect_each_orientation"]
    axes[1, 0].axvline(min_effect, linestyle="--", color="black", linewidth=1)
    axes[1, 0].set(xlabel="Minimum intended effect across orientations", ylabel="Variants", title="Robust effect distribution")
    axes[1, 1].axhline(min_effect, linestyle="--", color="black", linewidth=1)
    axes[1, 1].set(xlabel="Original candidate rank", ylabel="Minimum intended effect across orientations", title="Effects across candidate ranking")
    fig.tight_layout()
    fig.savefig(output / "perturbation_full_qc.png", dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
