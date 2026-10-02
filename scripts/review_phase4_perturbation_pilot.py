#!/usr/bin/env python3
"""Validate pilot sequence edits and summarize direct model effects."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    table = pd.read_csv(args.results, sep="\t")
    assert len(table) == 1000 and table.peak_id.nunique() == 100
    assert table.variant_id.nunique() == 1000
    assert np.isfinite(table.select_dtypes("number").to_numpy()).all()
    for row in table.itertuples():
        assert len(row.parent_sequence) == len(row.variant_sequence) == 500
        changed = [i for i, (a, b) in enumerate(zip(row.parent_sequence, row.variant_sequence)) if a != b]
        assert changed == [row.parent_position_0based]
        assert row.parent_sequence[changed[0]] == row.ref
        assert row.variant_sequence[changed[0]] == row.alt
        assert row.genomic_position_0based == row.start + changed[0]
        for field in ("parent", "variant"):
            assert hashlib.sha256(getattr(row, field + "_sequence").encode()).hexdigest() == getattr(row, field + "_sequence_sha256")

    # Project screening thresholds for this review; these do not prove biological activity.
    min_effect = 0.1
    max_profile_jsd = 0.05
    profile_ok = (table.forward_profile_jsd_vs_parent <= max_profile_jsd) & (table.reverse_profile_jsd_vs_parent <= max_profile_jsd)
    table["pilot_edit_screen_pass"] = (table.robust_effect_lower_bound >= min_effect) & profile_ok
    table["pilot_screen_min_effect_each_orientation"] = min_effect
    table["pilot_screen_max_profile_jsd_each_orientation"] = max_profile_jsd
    summary = {"status": "GO_FOR_SCALED_SINGLE_BASE_NOMINATION_AND_DIRECT_RESCORING", "parents": 100,
               "variants": 1000, "min_effect_each_orientation": min_effect,
               "max_profile_jsd_each_orientation": max_profile_jsd, "by_design_class": {}}
    for design_class, group in table.groupby("design_class"):
        selected = group[group.pilot_edit_screen_pass]
        summary["by_design_class"][design_class] = {
            "variants": len(group),
            "intended_effect_both_orientations": int(group.intended_effect_in_both_orientations.sum()),
            "screen_pass_variants": len(selected), "screen_pass_parents": selected.peak_id.nunique(),
            "median_delta_mean_logcount": float(group.delta_mean_logcount.median()),
            "median_robust_effect_lower_bound": float(group.robust_effect_lower_bound.median()),
            "nomination_direct_effect_spearman": float(group.nomination_count_hypothetical_delta.corr(group.delta_mean_logcount, method="spearman")),
            "forward_reverse_effect_pearson": float(group.delta_forward_logcount.corr(group.delta_reverse_complement_logcount)),
            "profile_jsd_failures": int((~profile_ok.loc[group.index]).sum()),
        }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "perturbation_pilot_review.tsv.gz", sep="\t", index=False, compression="gzip")
    (args.output_dir / "pilot_review_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for name, color in (("gain", "#4c78a8"), ("loss_control", "#e45756")):
        group = table[table.design_class == name]
        axes[0, 0].scatter(group.delta_forward_logcount, group.delta_reverse_complement_logcount, s=10, alpha=.5, color=color, label=name)
        axes[0, 1].scatter(group.nomination_count_hypothetical_delta, group.delta_mean_logcount, s=10, alpha=.5, color=color, label=name)
        axes[1, 0].hist(group.robust_effect_lower_bound, bins=40, alpha=.6, color=color, label=name)
        axes[1, 1].hist(group.mean_profile_jsd_vs_parent, bins=40, alpha=.6, color=color, label=name)
    axes[0, 0].axhline(0, color="black", linewidth=.8)
    axes[0, 0].axvline(0, color="black", linewidth=.8)
    axes[0, 0].set(xlabel="Forward count effect (log units)", ylabel="RC count effect (log units)", title="Agreement of directly scored effects")
    axes[0, 1].set(xlabel="Hypothetical attribution nomination delta", ylabel="Direct mean count effect (log units)", title="Attribution calibration")
    axes[1, 0].axvline(min_effect, linestyle="--", color="black", linewidth=1)
    axes[1, 0].set(xlabel="Minimum intended effect across orientations", ylabel="Variants", title="Robust predicted effect")
    axes[1, 1].set(xlabel="Mean profile JSD vs parent", ylabel="Variants", title="Profile change after substitution")
    axes[0, 0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.output_dir / "perturbation_pilot_qc.png", dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
