#!/usr/bin/env python3
"""Assemble and quality-control the Phase 3 donor pseudobulk count matrix."""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DONORS = ("LGS1", "LGS2", "LGS3", "LVG1")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(2 ** 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_peak_ids(path):
    rows = []
    with path.open() as handle:
        for line in handle:
            chrom, start, end, peak_id = line.rstrip("\n").split("\t")
            rows.append((peak_id, chrom, int(start), int(end)))
    return rows


def load_counts(path, expected_ids):
    values = np.empty(len(expected_ids), dtype=np.int64)
    with gzip.open(str(path), "rt") as handle:
        header = handle.readline().rstrip("\n")
        if header != "peak_id\tcount":
            raise ValueError("Unexpected count header in {}".format(path))
        count = 0
        for count, line in enumerate(handle, 1):
            peak_id, value = line.rstrip("\n").split("\t")
            if count > len(expected_ids) or peak_id != expected_ids[count - 1]:
                raise ValueError("Peak-order mismatch in {} at row {}".format(path, count))
            values[count - 1] = int(value)
    if count != len(expected_ids):
        raise ValueError("Expected {} peaks in {}; found {}".format(len(expected_ids), path, count))
    return values


def eta_squared(values, labels):
    values = np.asarray(values, dtype=float)
    grand = values.mean()
    total = np.square(values - grand).sum()
    if total == 0:
        return 0.0
    between = 0.0
    for label in sorted(set(labels)):
        subset = values[np.asarray(labels) == label]
        between += len(subset) * float((subset.mean() - grand) ** 2)
    return float(between / total)


def write_tsv(path, fieldnames, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def plot_sample_qc(rows, groups, output):
    figure, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    palette = plt.get_cmap("tab20")
    colors = {group: palette(index / max(1, len(groups) - 1)) for index, group in enumerate(groups)}
    x = np.arange(len(rows))
    axes[0].bar(x, [row["fragment_records"] / 1e6 for row in rows], color=[colors[row["analysis_group"]] for row in rows])
    axes[0].set_ylabel("Retained fragments (millions)")
    axes[1].bar(x, [row["frip"] for row in rows], color=[colors[row["analysis_group"]] for row in rows])
    axes[1].set_ylabel("FRiP")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([row["sample_id"] for row in rows], rotation=90, fontsize=7)
    for axis in axes:
        axis.grid(axis="y", alpha=0.25)
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_correlation(correlation, sample_ids, output):
    figure, axis = plt.subplots(figsize=(13, 11))
    image = axis.imshow(correlation, vmin=-1, vmax=1, cmap="coolwarm", interpolation="nearest")
    axis.set_xticks(range(len(sample_ids)))
    axis.set_yticks(range(len(sample_ids)))
    axis.set_xticklabels(sample_ids, rotation=90, fontsize=5)
    axis.set_yticklabels(sample_ids, fontsize=5)
    for boundary in range(4, len(sample_ids), 4):
        axis.axhline(boundary - 0.5, color="white", linewidth=0.25)
        axis.axvline(boundary - 0.5, color="white", linewidth=0.25)
    figure.colorbar(image, ax=axis, label="Pearson r (log2 CPM)")
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_pca(scores, groups, donors, explained, output):
    figure, axis = plt.subplots(figsize=(10, 8))
    palette = plt.get_cmap("tab20")
    group_order = list(dict.fromkeys(groups))
    colors = {group: palette(index / max(1, len(group_order) - 1)) for index, group in enumerate(group_order)}
    markers = {"LGS1": "o", "LGS2": "s", "LGS3": "^", "LVG1": "D"}
    for group in group_order:
        for donor in DONORS:
            selected = [index for index, value in enumerate(groups) if value == group and donors[index] == donor]
            axis.scatter(scores[selected, 0], scores[selected, 1], color=colors[group], marker=markers[donor], s=55, edgecolor="black", linewidth=0.25)
    for index, (x, y) in enumerate(scores[:, :2]):
        axis.annotate(donors[index], (x, y), xytext=(3, 2), textcoords="offset points", fontsize=5)
    handles = [plt.Line2D([], [], color=colors[group], marker="o", linestyle="", label=group) for group in group_order]
    axis.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
    axis.set_xlabel("PC1 ({:.1f}%)".format(100 * explained[0]))
    axis.set_ylabel("PC2 ({:.1f}%)".format(100 * explained[1]))
    axis.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count-root", type=Path, required=True)
    parser.add_argument("--peaks", type=Path, required=True)
    parser.add_argument("--cell-types", type=Path, required=True)
    parser.add_argument("--pseudobulk-qc", type=Path, required=True)
    parser.add_argument("--matrix-output-dir", type=Path, required=True)
    parser.add_argument("--report-output-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.matrix_output_dir.exists():
        raise FileExistsError("Refusing to overwrite: {}".format(args.matrix_output_dir))
    if args.report_output_dir.exists():
        raise FileExistsError("Refusing to overwrite: {}".format(args.report_output_dir))
    args.matrix_output_dir.mkdir(parents=True)
    args.report_output_dir.mkdir(parents=True)

    cell_types = read_tsv(args.cell_types)
    groups = [row["analysis_group"] for row in cell_types]
    sample_rows = [
        {"sample_id": "{}__{}".format(group, donor), "analysis_group": group, "donor": donor}
        for group in groups
        for donor in DONORS
    ]
    expected_depth = {
        (row["analysis_group"], row["donor"]): int(row["retained_fragment_records"])
        for row in read_tsv(args.pseudobulk_qc)
    }
    peaks = load_peak_ids(args.peaks)
    peak_ids = [row[0] for row in peaks]
    peak_sha = sha256(args.peaks)
    matrix = np.empty((len(peaks), len(sample_rows)), dtype=np.int64)
    qc_rows = []
    for column, sample in enumerate(sample_rows):
        sample_dir = args.count_root / sample["sample_id"]
        summary_path = sample_dir / "count_summary.json"
        count_path = sample_dir / "peak_counts.tsv.gz"
        checksum_path = sample_dir / "output_checksums.sha256"
        for path in (summary_path, count_path, checksum_path):
            if not path.is_file():
                raise FileNotFoundError(path)
        for line in checksum_path.read_text().splitlines():
            checksum, filename = line.split(None, 1)
            filename = filename.strip()
            if sha256(sample_dir / filename) != checksum:
                raise ValueError("Checksum failure: {}".format(sample_dir / filename))
        summary = json.loads(summary_path.read_text())
        if summary["sample_id"] != sample["sample_id"] or summary["peak_file_sha256"] != peak_sha:
            raise ValueError("Sample or peak-universe mismatch: {}".format(sample["sample_id"]))
        if summary["fragment_records"] != expected_depth[(sample["analysis_group"], sample["donor"])]:
            raise ValueError("Pseudobulk-depth mismatch: {}".format(sample["sample_id"]))
        matrix[:, column] = load_counts(count_path, peak_ids)
        if int(matrix[:, column].sum()) != summary["peak_overlap_assignments"]:
            raise ValueError("Count-sum mismatch: {}".format(sample["sample_id"]))
        qc_rows.append(dict(sample, **{key: summary[key] for key in (
            "fragment_records", "fragments_in_any_peak", "frip", "peak_overlap_assignments",
            "fragments_overlapping_multiple_peaks", "extra_peak_assignments", "detected_peaks",
            "zero_count_peaks", "maximum_peak_count")}))

    matrix_path = args.matrix_output_dir / "retina_consensus_peak_counts.tsv.gz"
    part_path = Path(str(matrix_path) + ".part")
    with gzip.open(str(part_path), "wt", compresslevel=1) as handle:
        handle.write("peak_id\tchrom\tstart\tend\t{}\n".format("\t".join(row["sample_id"] for row in sample_rows)))
        for row_index, peak in enumerate(peaks):
            handle.write("{}\t{}\t{}\t{}\t{}\n".format(
                peak[0], peak[1], peak[2], peak[3],
                "\t".join(str(value) for value in matrix[row_index, :]),
            ))
    part_path.replace(matrix_path)
    (args.matrix_output_dir / "matrix.sha256").write_text(
        "{}  {}\n".format(sha256(matrix_path), matrix_path.name)
    )

    library_sizes = matrix.sum(axis=0).astype(float)
    if np.any(library_sizes <= 0):
        raise ValueError("At least one sample has an empty peak-count library")
    cpm = matrix.astype(float) / library_sizes * 1e6
    exploratory_filter = np.sum(cpm >= 1.0, axis=1) >= len(DONORS)
    if int(np.sum(exploratory_filter)) < 2:
        raise ValueError("Fewer than two peaks pass the exploratory QC filter")
    log_cpm = np.log2(
        (matrix[exploratory_filter, :].astype(float) + 0.5)
        / (library_sizes + 1.0)
        * 1e6
    )
    correlation = np.corrcoef(log_cpm.T)
    variance = np.var(log_cpm, axis=1)
    variable = np.flatnonzero(variance > 0)
    variable = variable[np.argsort(variance[variable])[-min(20000, len(variable)):]]
    pca_input = log_cpm[variable, :].T
    pca_input -= pca_input.mean(axis=0, keepdims=True)
    u, singular, _ = np.linalg.svd(pca_input, full_matrices=False)
    scores = u * singular
    explained = np.square(singular) / np.square(singular).sum()
    sample_groups = [row["analysis_group"] for row in sample_rows]
    sample_donors = [row["donor"] for row in sample_rows]

    for index, row in enumerate(qc_rows):
        row["peak_count_library_size"] = int(library_sizes[index])
        row["detected_peaks_cpm_ge_1"] = int(np.sum(cpm[:, index] >= 1))
    qc_fields = list(qc_rows[0])
    write_tsv(args.report_output_dir / "sample_qc.tsv", qc_fields, qc_rows)

    correlation_rows = []
    sample_ids = [row["sample_id"] for row in sample_rows]
    for index, sample_id in enumerate(sample_ids):
        correlation_rows.append(dict(sample_id=sample_id, **{
            other: float(correlation[index, column]) for column, other in enumerate(sample_ids)
        }))
    write_tsv(args.report_output_dir / "sample_correlation_logcpm.tsv", ["sample_id"] + sample_ids, correlation_rows)

    pca_rows = []
    for index, sample in enumerate(sample_rows):
        row = dict(sample)
        for component in range(min(10, scores.shape[1])):
            row["PC{}".format(component + 1)] = float(scores[index, component])
        pca_rows.append(row)
    write_tsv(args.report_output_dir / "pca_scores.tsv", list(pca_rows[0]), pca_rows)

    within = {}
    for group_index, group in enumerate(groups):
        indexes = list(range(group_index * 4, group_index * 4 + 4))
        values = [correlation[a, b] for offset, a in enumerate(indexes) for b in indexes[offset + 1:]]
        within[group] = {
            "minimum": float(np.min(values)),
            "median": float(np.median(values)),
            "maximum": float(np.max(values)),
        }
    associations = []
    for component in range(min(10, scores.shape[1])):
        associations.append({
            "component": component + 1,
            "explained_variance_fraction": float(explained[component]),
            "cell_type_eta_squared": eta_squared(scores[:, component], sample_groups),
            "donor_eta_squared": eta_squared(scores[:, component], sample_donors),
        })
    result = {
        "status": "COUNT_MATRIX_QC_READY_FOR_SCIENTIFIC_REVIEW",
        "differential_accessibility_status": "BLOCKED_PENDING_COUNT_MATRIX_QC_REVIEW",
        "samples": len(sample_rows),
        "consensus_peaks": len(peaks),
        "matrix_sha256": sha256(matrix_path),
        "total_fragment_records": int(sum(row["fragment_records"] for row in qc_rows)),
        "total_fragments_in_any_peak": int(sum(row["fragments_in_any_peak"] for row in qc_rows)),
        "frip_range": [float(min(row["frip"] for row in qc_rows)), float(max(row["frip"] for row in qc_rows))],
        "detected_peak_range": [int(min(row["detected_peaks"] for row in qc_rows)), int(max(row["detected_peaks"] for row in qc_rows))],
        "multiple_peak_fragment_fraction": float(
            sum(row["fragments_overlapping_multiple_peaks"] for row in qc_rows)
            / sum(row["fragment_records"] for row in qc_rows)
        ),
        "exploratory_qc_peaks": int(np.sum(exploratory_filter)),
        "exploratory_qc_filter": "CPM >= 1 in at least 4 of 52 donor pseudobulks",
        "within_cell_type_donor_correlations": within,
        "pca_associations": associations,
        "normalization_scope": "Exploratory log2 CPM only; edgeR normalization factors are not estimated here.",
        "review_required": [
            "FRiP and detected-peak outliers",
            "within-cell-type donor concordance",
            "cell-type separation versus donor association in PCA",
            "low-cell astrocyte, microglia, and RGC strata",
            "rod depth dominance",
        ],
    }
    (args.report_output_dir / "count_matrix_qc.json").write_text(json.dumps(result, indent=2) + "\n")
    plot_sample_qc(qc_rows, groups, args.report_output_dir / "sample_depth_frip.png")
    plot_correlation(correlation, sample_ids, args.report_output_dir / "sample_correlation_logcpm.png")
    plot_pca(scores, sample_groups, sample_donors, explained, args.report_output_dir / "pca_logcpm.png")
    with (args.report_output_dir / "output_checksums.sha256").open("w") as handle:
        for path in sorted(args.report_output_dir.iterdir()):
            if path.name != "output_checksums.sha256":
                handle.write("{}  {}\n".format(sha256(path), path.name))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
