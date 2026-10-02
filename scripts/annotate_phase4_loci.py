#!/usr/bin/env python3
"""Annotate natural parent loci and proposed edit positions with hg38 tracks."""
import argparse
from bisect import bisect_left, bisect_right
import gzip
import json
from pathlib import Path
import re

from download_phase4_annotation_assets import digest


def union_bp(intervals):
    total = 0
    end = None
    for start, stop in sorted(intervals):
        if end is None or start >= end:
            total += stop - start
        elif stop > end:
            total += stop - end
        end = stop if end is None else max(end, stop)
    return total


def repeat_features(hits, start, end):
    clipped = [(max(start, hit[0]), min(end, hit[1])) for hit in hits
               if hit[0] < end and hit[1] > start]
    overlapping = [hit for hit in hits if hit[0] < end and hit[1] > start]
    covered = union_bp(clipped)
    return {
        "repeat_overlap_bp": covered, "repeat_overlap_fraction": covered / (end - start),
        "repeat_hit_count": len(overlapping),
        "repeat_names": ";".join(sorted({hit[2] for hit in overlapping})),
        "repeat_classes": ";".join(sorted({hit[3] for hit in overlapping})),
        "repeat_families": ";".join(sorted({hit[4] for hit in overlapping})),
    }


def mappability_features(values):
    import numpy as np
    values = np.asarray(values, dtype=float)
    if not len(values) or np.isinf(values).any():
        raise ValueError("Invalid mappability values")
    omitted = np.isnan(values)
    values = np.where(omitted, 0.0, values)
    if np.any(values < 0) or np.any(values > 1 + 1e-6):
        raise ValueError("Mappability values outside [0,1]")
    return {
        "mean": float(values.mean()), "min": float(values.min()),
        "zero_fraction": float((values == 0).mean()),
        "fraction_ge_0_9": float((values >= .9).mean()),
        "omitted_fraction": float(omitted.mean()),
    }, values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    import pandas as pd
    import pyBigWig

    config = json.loads(args.config.read_text())
    root = Path(config["assets_directory"])
    manifest = json.loads((root / "asset_snapshot.json").read_text())
    for asset in config["assets"]:
        record = manifest["assets"][asset["name"]]
        if record["url"] != asset["url"] or digest(root / asset["name"]) != record["sha256"]:
            raise ValueError("Annotation asset checksum or URL mismatch")
    sql = (root / "rmsk.sql").read_text()
    columns = re.findall(r"^\s+`(\w+)`", sql, flags=re.MULTILINE)
    expected = ["bin", "swScore", "milliDiv", "milliDel", "milliIns", "genoName",
                "genoStart", "genoEnd", "genoLeft", "strand", "repName", "repClass",
                "repFamily", "repStart", "repEnd", "repLeft", "id"]
    if columns != expected:
        raise ValueError("Unexpected UCSC rmsk schema")

    parents = pd.read_csv(config["parents"], sep="\t")
    edits = pd.read_csv(config["edits"], sep="\t")
    if len(parents) != 1000 or parents.peak_id.nunique() != 1000:
        raise ValueError("Expected all 1,000 original natural parents")
    if edits.variant_id.duplicated().any():
        raise ValueError("Duplicate proposed variant identifiers")
    lengths = {line.split()[0]: int(line.split()[1]) for line in Path(config["reference_fai"]).read_text().splitlines()}
    if any(int(row.end) - int(row.start) != 500 for row in parents.itertuples()):
        raise ValueError("Unexpected parent length")

    # Retain only repeats that intersect at least one selected parent.
    by_chrom = {}
    for chrom, group in parents.groupby("chrom"):
        group = group.sort_values("start")
        by_chrom[chrom] = (group.start.to_list(), list(group.itertuples()))
    hits = {peak: [] for peak in parents.peak_id}
    with gzip.open(root / "rmsk.txt.gz", "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(expected):
                raise ValueError("Malformed rmsk row")
            chrom = fields[5]
            if chrom not in by_chrom:
                continue
            start, end = int(fields[6]), int(fields[7])
            positions, rows = by_chrom[chrom]
            lower = bisect_right(positions, start - 500)
            upper = bisect_left(positions, end)
            for row in rows[lower:upper]:
                hits[row.peak_id].append((start, end, fields[10], fields[11], fields[12]))

    output = parents.copy()
    features = [repeat_features(hits[row.peak_id], int(row.start), int(row.end))
                for row in parents.itertuples()]
    output = pd.concat([output.reset_index(drop=True), pd.DataFrame(features)], axis=1)
    base_values = {}
    for k in config["mappability"]["read_lengths_bp"]:
        with pyBigWig.open(str(root / "k{}.Umap.MultiTrackMappability.bw".format(k))) as bw:
            chroms = bw.chroms()
            for chrom in parents.chrom.unique():
                if chrom not in chroms or chroms[chrom] != lengths[chrom]:
                    raise ValueError("Mappability/reference chromosome mismatch: {}".format(chrom))
            stats = []
            for row in parents.itertuples():
                values = bw.values(row.chrom, int(row.start), int(row.end), numpy=True)
                summary, values = mappability_features(values)
                stats.append({"umap_k{}_parent_{}".format(k, name): value for name, value in summary.items()})
                base_values[(row.peak_id, k)] = values
        output = pd.concat([output, pd.DataFrame(stats)], axis=1)

    parent_index = parents.set_index("peak_id")
    edit_features = []
    for row in edits.itertuples():
        parent = parent_index.loc[row.peak_id]
        position = int(row.parent_position_0based)
        if (row.chrom != parent.chrom or row.start != parent.start or row.end != parent.end
                or row.parent_sequence != parent.parent_sequence
                or not 0 <= position < 500 or row.genomic_position_0based != row.start + position):
            raise ValueError("Edit/parent locus mismatch")
        features = repeat_features(hits[row.peak_id], row.genomic_position_0based, row.genomic_position_0based + 1)
        features = {"edit_base_" + key: value for key, value in features.items()}
        for k in config["mappability"]["read_lengths_bp"]:
            features["edit_base_reference_umap_k{}".format(k)] = float(base_values[(row.peak_id, k)][position])
        features["variant_mappability_prediction_available"] = False
        edit_features.append(features)
    annotated_edits = pd.concat([edits.reset_index(drop=True), pd.DataFrame(edit_features)], axis=1)
    annotation_columns = [column for column in output.columns if column not in parents.columns]
    annotated_edits = annotated_edits.merge(
        output[["peak_id"] + annotation_columns], on="peak_id", how="left", validate="many_to_one"
    )

    args.output_dir.mkdir(parents=True, exist_ok=False)
    output.to_csv(args.output_dir / "annotated_natural_parents.tsv.gz", sep="\t", index=False, compression="gzip")
    annotated_edits.to_csv(args.output_dir / "annotated_best_single_edits.tsv.gz", sep="\t", index=False, compression="gzip")
    repeat_rows = [{"peak_id": peak, "repeat_start": hit[0], "repeat_end": hit[1],
                    "repeat_name": hit[2], "repeat_class": hit[3], "repeat_family": hit[4]}
                   for peak, values in hits.items() for hit in values]
    pd.DataFrame(repeat_rows, columns=["peak_id", "repeat_start", "repeat_end", "repeat_name", "repeat_class", "repeat_family"]).to_csv(
        args.output_dir / "parent_repeat_hits.tsv.gz", sep="\t", index=False, compression="gzip")
    summary = {
        "status": "LOCUS_ANNOTATIONS_READY_FOR_REVIEW", "parents": len(output), "proposed_edits": len(annotated_edits),
        "parents_overlapping_repeats": int((output.repeat_overlap_bp > 0).sum()),
        "edits_in_repeats": int((annotated_edits.edit_base_repeat_overlap_bp > 0).sum()),
        "repeat_overlap_fraction_quantiles": {str(q): float(output.repeat_overlap_fraction.quantile(q)) for q in (0,.25,.5,.75,.95,1)},
        "reference_locus_mappability": {}, "automatic_exclusions": 0,
        "next_gate": config["next_gate"],
    }
    for k in config["mappability"]["read_lengths_bp"]:
        series = output["umap_k{}_parent_mean".format(k)]
        summary["reference_locus_mappability"][str(k)] = {
            "mean_quantiles": {str(q): float(series.quantile(q)) for q in (0,.05,.25,.5,.75,.95,1)},
            "parents_mean_lt_0_8_descriptive": int((series < .8).sum()),
        }
    (args.output_dir / "annotation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    manifest["analysis_inputs"] = {name: {"path": config[name], "sha256": digest(Path(config[name]))}
                                   for name in ("parents", "edits", "reference_fai")}
    (args.output_dir / "annotation_provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
