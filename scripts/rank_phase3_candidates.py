#!/usr/bin/env python3
"""Rank Muller-glia-accessible regions after the differential-model QC gate."""

import argparse
from array import array
import bisect
import csv
import gzip
import json
import math
from pathlib import Path


DONORS = ("LGS1", "LGS2", "LGS3", "LVG1")


def open_text(path, mode="rt"):
    path = Path(path)
    return gzip.open(path, mode, newline="") if path.suffix == ".gz" else path.open(mode, newline="")


def header_indexes(header):
    """Return the last index for each name; edgeR outputs repeat coordinate names."""
    return {name: index for index, name in enumerate(header)}


def read_primary(path):
    peaks = []
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        index = header_indexes(header)
        required = ["peak_id", "chrom", "start", "end", "logFC", "logCPM", "FDR",
                    "donors_positive", "minimum_donor_logFC"]
        for name in required:
            if name not in index:
                raise ValueError(f"Missing column {name} in {path}")
        for row in reader:
            peaks.append({
                "peak_id": row[index["peak_id"]],
                "chrom": row[index["chrom"]],
                "start": int(row[index["start"]]),
                "end": int(row[index["end"]]),
                "primary_logFC": float(row[index["logFC"]]),
                "primary_logCPM": float(row[index["logCPM"]]),
                "primary_FDR": float(row[index["FDR"]]),
                "primary_donors_positive": int(row[index["donors_positive"]]),
                "primary_minimum_donor_logFC": float(row[index["minimum_donor_logFC"]]),
            })
    if not peaks or len({peak["peak_id"] for peak in peaks}) != len(peaks):
        raise ValueError("Primary result is empty or contains duplicate peak IDs")
    return peaks


def read_ordered_metrics(path, expected_ids, fields):
    values = {field: array("d") for field in fields}
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        index = header_indexes(header)
        for name in ("peak_id", *fields):
            if name not in index:
                raise ValueError(f"Missing column {name} in {path}")
        count = 0
        for count, row in enumerate(reader, start=1):
            position = count - 1
            if position >= len(expected_ids) or row[index["peak_id"]] != expected_ids[position]:
                raise ValueError(f"Peak order mismatch in {path} at row {count}")
            for field in fields:
                values[field].append(float(row[index[field]]))
    if count != len(expected_ids):
        raise ValueError(f"Expected {len(expected_ids)} rows in {path}, observed {count}")
    return values


def load_effective_libraries(path):
    with open_text(path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    effective = {row["sample_id"]: float(row["effective_library_size"]) for row in rows}
    if len(effective) != 52 or any(value <= 0 for value in effective.values()):
        raise ValueError("Expected 52 positive effective library sizes")
    return effective


def calculate_accessibility(matrix_path, peaks, effective):
    peak_index = {peak["peak_id"]: index for index, peak in enumerate(peaks)}
    target = array("d", [math.nan]) * len(peaks)
    maximum = array("d", [math.nan]) * len(peaks)
    maximum_group = [None] * len(peaks)
    with open_text(matrix_path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[:4] != ["peak_id", "chrom", "start", "end"]:
            raise ValueError("Unexpected count-matrix coordinate columns")
        samples = header[4:]
        if set(samples) != set(effective):
            raise ValueError("Count-matrix samples do not match normalization table")
        groups = {}
        for column, sample in enumerate(samples, start=4):
            group = sample.split("__", 1)[0]
            groups.setdefault(group, []).append((column, effective[sample]))
        if "MG" not in groups or len(groups) != 13 or any(len(columns) != 4 for columns in groups.values()):
            raise ValueError("Expected four donors for each of 13 cell types")
        found = 0
        for row in reader:
            position = peak_index.get(row[0])
            if position is None:
                continue
            means = {}
            for group, columns in groups.items():
                means[group] = sum(float(row[column]) / library * 1e6 for column, library in columns) / len(columns)
            off_target, off_value = max(
                ((group, value) for group, value in means.items() if group != "MG"),
                key=lambda item: (item[1], item[0]),
            )
            target[position] = means["MG"]
            maximum[position] = off_value
            maximum_group[position] = off_target
            found += 1
    if found != len(peaks) or any(math.isnan(value) for value in target):
        raise ValueError(f"Accessibility found for {found} of {len(peaks)} tested peaks")
    return target, maximum, maximum_group


def load_refgene(path):
    records = {}
    tss = {}
    canonical = {f"chr{i}" for i in range(1, 23)} | {"chrX"}
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 13 or row[2] not in canonical:
                continue
            chrom, strand = row[2], row[3]
            tx_start, tx_end = int(row[4]), int(row[5])
            exon_starts = tuple(int(value) for value in row[9].rstrip(",").split(",") if value)
            exon_ends = tuple(int(value) for value in row[10].rstrip(",").split(",") if value)
            gene = row[12] or row[1]
            gene_tss = tx_start if strand == "+" else tx_end
            record = (tx_start, tx_end, gene, gene_tss, tuple(zip(exon_starts, exon_ends)))
            records.setdefault(chrom, []).append(record)
            tss.setdefault(chrom, []).append((gene_tss, gene))
    for chrom in records:
        records[chrom].sort()
        tss[chrom].sort()
    return records, tss


def overlaps(start, end, other_start, other_end):
    return start < other_end and other_start < end


def annotate_peaks(peaks, refgene_path, promoter_window=2000):
    transcripts, tss = load_refgene(refgene_path)
    annotations = [None] * len(peaks)
    by_chrom = {}
    for index, peak in enumerate(peaks):
        by_chrom.setdefault(peak["chrom"], []).append(index)
    for chrom, indices in by_chrom.items():
        tx = transcripts.get(chrom, [])
        tss_records = tss.get(chrom, [])
        tss_positions = [item[0] for item in tss_records]
        cursor = 0
        active = []
        for index in sorted(indices, key=lambda value: peaks[value]["start"]):
            peak = peaks[index]
            start, end = peak["start"], peak["end"]
            while cursor < len(tx) and tx[cursor][0] < end:
                active.append(tx[cursor]); cursor += 1
            active = [record for record in active if record[1] > start]
            midpoint = (start + end) // 2
            position = bisect.bisect_left(tss_positions, midpoint)
            nearby = []
            for candidate in (position - 1, position):
                if 0 <= candidate < len(tss_records):
                    gene_tss, gene = tss_records[candidate]
                    nearby.append((abs(midpoint - gene_tss), gene, midpoint - gene_tss))
            nearest_distance, nearest_gene, signed_distance = min(nearby) if nearby else (None, "", None)
            promoter_left = bisect.bisect_left(tss_positions, start - promoter_window)
            promoter_right = bisect.bisect_right(tss_positions, end + promoter_window)
            promoter_genes = sorted({gene for _, gene in tss_records[promoter_left:promoter_right]})
            exonic = [record for record in active if any(overlaps(start, end, a, b) for a, b in record[4])]
            if promoter_genes:
                genomic_class = "promoter"
                genes = promoter_genes
            elif exonic:
                genomic_class = "exonic"
                genes = sorted({record[2] for record in exonic})
            elif active:
                genomic_class = "intronic"
                genes = sorted({record[2] for record in active})
            else:
                genomic_class = "distal_intergenic"
                genes = []
            annotations[index] = (genomic_class, ",".join(genes), nearest_gene, signed_distance)
    if any(value is None for value in annotations):
        raise ValueError("Some peaks were not annotated")
    return annotations


def percentiles(values):
    if not values:
        return []
    order = sorted(range(len(values)), key=lambda index: (values[index], index))
    result = [0.0] * len(values)
    denominator = max(1, len(values) - 1)
    begin = 0
    while begin < len(order):
        end = begin + 1
        while end < len(order) and values[order[end]] == values[order[begin]]:
            end += 1
        rank = ((begin + end - 1) / 2) / denominator
        for ordered_index in order[begin:end]:
            result[ordered_index] = rank
        begin = end
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--count-matrix", required=True, type=Path)
    parser.add_argument("--normalization", required=True, type=Path)
    parser.add_argument("--refgene", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    target_count = int(config["target_count"])
    contrasts_dir = args.results_dir / "contrasts"
    lodo_dir = args.results_dir / "leave_one_donor_out"
    peaks = read_primary(contrasts_dir / "MG_vs_rest.tsv.gz")
    peak_ids = [peak["peak_id"] for peak in peaks]

    contrast_names = sorted(path.name[len("MG_vs_"):-len(".tsv.gz")]
                            for path in contrasts_dir.glob("MG_vs_*.tsv.gz") if path.name != "MG_vs_rest.tsv.gz")
    if len(contrast_names) != 12:
        raise ValueError(f"Expected 12 pairwise contrasts, observed {len(contrast_names)}")
    pairwise = {}
    for name in contrast_names:
        pairwise[name] = read_ordered_metrics(
            contrasts_dir / f"MG_vs_{name}.tsv.gz", peak_ids,
            ("logFC", "FDR", "donors_positive", "minimum_donor_logFC"),
        )

    lodo = {}
    for comparator in ("Astrocyte", "Microglia"):
        for donor in DONORS:
            key = f"{comparator}_without_{donor}"
            lodo[key] = read_ordered_metrics(
                lodo_dir / f"MG_vs_{comparator}_without_{donor}.tsv.gz",
                peak_ids, ("logFC", "FDR"),
            )

    effective = load_effective_libraries(args.normalization)
    mg_cpm, max_off_cpm, max_off_group = calculate_accessibility(args.count_matrix, peaks, effective)
    annotations = annotate_peaks(peaks, args.refgene, config["annotation"]["promoter_window_bp"])

    eligibility = config["primary_eligibility"]
    eligible_indices = []
    raw_components = {name: [] for name in config["score_weights"]}
    derived = {}
    for index, peak in enumerate(peaks):
        pair_effects = [pairwise[name]["logFC"][index] for name in contrast_names]
        lodo_effects = [values["logFC"][index] for values in lodo.values()]
        lodo_fdrs = [values["FDR"][index] for values in lodo.values()]
        specificity = math.log2((mg_cpm[index] + 0.25) / (max_off_cpm[index] + 0.25))
        values = {
            "minimum_pairwise_log2fc_percentile": min(pair_effects),
            "observed_specificity_percentile": specificity,
            "primary_log2fc_percentile": peak["primary_logFC"],
            "primary_significance_percentile": -math.log10(max(peak["primary_FDR"], 1e-300)),
            "muller_accessibility_percentile": math.log2(mg_cpm[index] + 0.25),
            "primary_donor_floor_percentile": peak["primary_minimum_donor_logFC"],
            "glial_lodo_effect_floor_percentile": min(lodo_effects),
            "glial_lodo_significance_percentile": -math.log10(max(max(lodo_fdrs), 1e-300)),
        }
        derived[index] = (min(pair_effects), specificity, min(lodo_effects), max(lodo_fdrs))
        is_eligible = (
            peak["primary_FDR"] < eligibility["fdr_lt"]
            and peak["primary_logFC"] >= eligibility["log2fc_gte"]
            and peak["primary_donors_positive"] == eligibility["donors_positive"]
            and pairwise["Astrocyte"]["logFC"][index] > 0
            and pairwise["Microglia"]["logFC"][index] > 0
            and all(value > 0 for value in lodo_effects)
        )
        if is_eligible:
            eligible_indices.append(index)
            for name, value in values.items():
                raw_components[name].append(value)

    if len(eligible_indices) < target_count:
        raise ValueError(f"Only {len(eligible_indices)} regions pass eligibility; need {target_count}")
    component_percentiles = {name: percentiles(values) for name, values in raw_components.items()}
    scores = {}
    for eligible_position, peak_index in enumerate(eligible_indices):
        scores[peak_index] = sum(
            config["score_weights"][name] * component_percentiles[name][eligible_position]
            for name in config["score_weights"]
        )
    ranked = sorted(eligible_indices, key=lambda index: (-scores[index], peaks[index]["peak_id"]))
    ranks = {index: rank for rank, index in enumerate(ranked, start=1)}
    selected = set(ranked[:target_count])

    args.output_dir.mkdir(parents=True, exist_ok=False)
    pair_columns = [item for name in contrast_names for item in (f"{name}_logFC", f"{name}_FDR")]
    lodo_columns = [item for key in lodo for item in (f"{key}_logFC", f"{key}_FDR")]
    columns = [
        "candidate_rank", "selected_top_target", "eligible", "ranking_score",
        "peak_id", "chrom", "start", "end", "genomic_class", "overlapping_genes",
        "nearest_gene", "signed_distance_to_nearest_tss", "primary_logFC", "primary_FDR",
        "primary_logCPM", "primary_donors_positive", "primary_minimum_donor_logFC",
        "minimum_pairwise_logFC", "muller_mean_cpm", "maximum_offtarget_mean_cpm",
        "maximum_offtarget_cell_type", "observed_specificity_log2_ratio",
        "minimum_glial_lodo_logFC", "maximum_glial_lodo_FDR", *pair_columns, *lodo_columns,
    ]

    def output_row(index):
        peak = peaks[index]
        genomic_class, genes, nearest_gene, distance = annotations[index]
        minimum_pairwise, specificity, minimum_lodo, maximum_lodo_fdr = derived[index]
        row = [
            ranks.get(index, ""), index in selected, index in scores,
            f"{scores[index]:.10g}" if index in scores else "", peak["peak_id"], peak["chrom"],
            peak["start"], peak["end"], genomic_class, genes, nearest_gene,
            "" if distance is None else distance, peak["primary_logFC"], peak["primary_FDR"],
            peak["primary_logCPM"], peak["primary_donors_positive"], peak["primary_minimum_donor_logFC"],
            minimum_pairwise, mg_cpm[index], max_off_cpm[index], max_off_group[index], specificity,
            minimum_lodo, maximum_lodo_fdr,
        ]
        for name in contrast_names:
            row.extend((pairwise[name]["logFC"][index], pairwise[name]["FDR"][index]))
        for key in lodo:
            row.extend((lodo[key]["logFC"][index], lodo[key]["FDR"][index]))
        return row

    all_path = args.output_dir / "all_scored_regions.tsv.gz"
    with open_text(all_path, "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for index in range(len(peaks)):
            writer.writerow(output_row(index))
    top_path = args.output_dir / "top_1000_candidates.tsv.gz"
    with open_text(top_path, "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for index in ranked[:target_count]:
            writer.writerow(output_row(index))
    with (args.output_dir / "top_1000_candidates.bed").open("w") as handle:
        for index in ranked[:target_count]:
            peak = peaks[index]
            handle.write(f'{peak["chrom"]}\t{peak["start"]}\t{peak["end"]}\t{peak["peak_id"]}\t{ranks[index]}\n')

    class_counts = {}
    for index in ranked[:target_count]:
        genomic_class = annotations[index][0]
        class_counts[genomic_class] = class_counts.get(genomic_class, 0) + 1
    summary = {
        "status": "CANDIDATE_RANKING_READY_FOR_REVIEW",
        "tested_regions": len(peaks),
        "eligible_regions": len(eligible_indices),
        "selected_regions": target_count,
        "top_set_genomic_classes": class_counts,
        "config": str(args.config),
        "outputs": {"all_scored": str(all_path), "top_candidates": str(top_path)},
        "interpretation": config["interpretation"],
    }
    (args.output_dir / "candidate_ranking_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
