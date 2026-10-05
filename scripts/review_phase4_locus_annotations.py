#!/usr/bin/env python3
"""Verify natural-locus annotations and produce descriptive review flags."""
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
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--original-parents", required=True, type=Path)
    parser.add_argument("--original-edits", required=True, type=Path)
    args = parser.parse_args()
    root = args.results_dir
    for line in (root / "output_checksums.sha256").read_text().splitlines():
        expected, name = line.split()
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected
    parents = pd.read_csv(root / "annotated_natural_parents.tsv.gz", sep="\t")
    edits = pd.read_csv(root / "annotated_best_single_edits.tsv.gz", sep="\t")
    original_parents = pd.read_csv(args.original_parents, sep="\t")
    original_edits = pd.read_csv(args.original_edits, sep="\t")
    assert len(parents) == parents.peak_id.nunique() == 1000
    assert len(edits) == edits.variant_id.nunique() == 1884
    for annotated, original, key in ((parents, original_parents, "peak_id"), (edits, original_edits, "variant_id")):
        annotated = annotated.set_index(key).loc[original[key]].reset_index()
        pd.testing.assert_frame_equal(annotated[original.columns], original.reset_index(drop=True))
    assert ((parents.repeat_overlap_bp >= 0) & (parents.repeat_overlap_bp <= 500)).all()
    assert np.allclose(parents.repeat_overlap_fraction, parents.repeat_overlap_bp / 500)
    hits = pd.read_csv(root / "parent_repeat_hits.tsv.gz", sep="\t")
    # Independent per-base reconstruction of the union for each natural parent.
    by_peak = {peak: group for peak, group in hits.groupby("peak_id")}
    for row in parents.itertuples():
        covered = np.zeros(500, dtype=bool)
        group = by_peak.get(row.peak_id)
        if group is not None:
            for hit in group.itertuples():
                left, right = max(row.start, hit.repeat_start), min(row.end, hit.repeat_end)
                assert left < right
                covered[left - row.start:right - row.start] = True
        assert int(covered.sum()) == row.repeat_overlap_bp
    for table in (parents, edits):
        for k in (50, 100):
            for metric in ("mean", "min", "zero_fraction", "fraction_ge_0_9", "omitted_fraction"):
                values = table["umap_k{}_parent_{}".format(k, metric)]
                assert np.isfinite(values).all() and values.between(0, 1 + 1e-6).all()
        table["reference_umap_k50_mean_lt_0_8_review_flag"] = table.umap_k50_parent_mean < .8
        table["reference_umap_k100_mean_lt_0_8_review_flag"] = table.umap_k100_parent_mean < .8
        table["repeat_overlap_ge_50pct_review_flag"] = table.repeat_overlap_fraction >= .5
    parent_index = parents.set_index("peak_id")
    for row in edits.itertuples():
        assert row.edit_base_repeat_overlap_bp in (0, 1)
        group = by_peak.get(row.peak_id)
        expected = False if group is None else bool(((group.repeat_start <= row.genomic_position_0based) & (group.repeat_end > row.genomic_position_0based)).any())
        assert bool(row.edit_base_repeat_overlap_bp) == expected
        assert row.repeat_overlap_fraction == parent_index.loc[row.peak_id, "repeat_overlap_fraction"]
        assert not row.variant_mappability_prediction_available
    parents.to_csv(root / "parent_annotation_review.tsv.gz", sep="\t", index=False, compression="gzip")
    edits.to_csv(root / "edit_annotation_review.tsv.gz", sep="\t", index=False, compression="gzip")
    low_map = parents[(parents.umap_k50_parent_mean < .8) | (parents.umap_k100_parent_mean < .8)]
    low_map.to_csv(root / "parents_with_low_reference_mappability.tsv.gz", sep="\t", index=False, compression="gzip")
    summary = {
        "status": "PASS_ANNOTATION_INTEGRITY_LIBRARY_REVIEW_PENDING",
        "parents": len(parents), "proposed_edits": len(edits),
        "parents_repeat_overlap_any": int((parents.repeat_overlap_bp > 0).sum()),
        "parents_repeat_overlap_ge_50pct": int((parents.repeat_overlap_fraction >= .5).sum()),
        "parents_repeat_overlap_ge_90pct": int((parents.repeat_overlap_fraction >= .9).sum()),
        "parents_umap_k50_mean_lt_0_8": int((parents.umap_k50_parent_mean < .8).sum()),
        "parents_umap_k100_mean_lt_0_8": int((parents.umap_k100_parent_mean < .8).sum()),
        "proposed_edits_in_repeats": int((edits.edit_base_repeat_overlap_bp > 0).sum()),
        "proposed_edit_parents_with_umap_k50_mean_lt_0_8": int(edits[edits.umap_k50_parent_mean < .8].peak_id.nunique()),
        "automatic_exclusions": 0, "by_edit_class": {},
        "next_gate": "ASSEMBLE_PARENT_VARIANT_LIBRARY_REVIEW_TABLE_WITH_MOTIF_EVIDENCE",
    }
    for kind, group in edits.groupby("design_class"):
        summary["by_edit_class"][kind] = {
            "edits": len(group), "edits_in_repeats": int((group.edit_base_repeat_overlap_bp > 0).sum()),
            "edit_base_reference_umap_k50_lt_0_8": int((group.edit_base_reference_umap_k50 < .8).sum()),
            "edit_base_reference_umap_k100_lt_0_8": int((group.edit_base_reference_umap_k100 < .8).sum()),
        }
    (root / "annotation_review_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    axes[0].hist(parents.repeat_overlap_fraction, bins=25, color="#4c78a8", edgecolor="white")
    axes[0].set(xlabel="Parent fraction overlapping repeats", ylabel="Parents", title="Repeat overlap across 1,000 parents")
    axes[1].scatter(parents.umap_k50_parent_mean, parents.umap_k100_parent_mean, c=parents.repeat_overlap_fraction, cmap="viridis", s=15, alpha=.65)
    axes[1].axvline(.8, linestyle="--", linewidth=1, color="black")
    axes[1].axhline(.8, linestyle="--", linewidth=1, color="black")
    axes[1].set(xlabel="50-bp mean Umap", ylabel="100-bp mean Umap", title="Natural reference mappability")
    scatter = axes[2].scatter(parents.repeat_overlap_fraction, parents.umap_k50_parent_mean, c=parents.repeat_overlap_fraction, cmap="viridis", s=15, alpha=.65)
    axes[2].set(xlabel="Parent repeat-overlap fraction", ylabel="50-bp mean Umap", title="Repeats and mapping context")
    fig.colorbar(scatter, ax=axes[2], label="Repeat-overlap fraction")
    fig.tight_layout()
    fig.savefig(root / "locus_annotation_qc.png", dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
