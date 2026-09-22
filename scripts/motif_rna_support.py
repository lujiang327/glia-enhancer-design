#!/usr/bin/env python3
"""Summarize paired-RNA support for TF families recovered by TF-MoDISco."""
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path


GENES = [
    "SOX9", "SOX10", "NFIA", "NFIB", "NFIC", "NFIX", "LHX2",
    "RFX1", "RFX2", "RFX3", "FOS", "FOSB", "JUN", "JUNB", "JUND",
    "CTCF", "RORA", "NFIL3", "DBP", "CRX", "OTX2", "PAX6", "NKX6-1",
    "NOTO", "VSX1", "MITF", "TFEB", "NRF1", "NFYA", "NFYB", "NFYC",
    "SP1", "SP3", "MEF2A", "MEF2B", "MEF2D", "GABPA", "IRF1", "IRF2",
    "REST", "YY1",
]


def run(root=Path(".")):
    raw = root / "data/raw"
    qc = root / "reports/qc"
    output = root / "reports/training/chrombpnet_neg1_motif_rna_support.tsv"
    samples = {
        row["rna_prefix"]: row
        for row in csv.DictReader(open(root / "config/samples.tsv"), delimiter="\t")
    }
    audit = {
        row["bare_barcode"]: row
        for row in csv.DictReader(open(qc / "barcode_audit.tsv"), delimiter="\t")
    }
    with gzip.open(raw / "GSM5866081_barcodes.tsv.gz", "rt") as handle:
        cells = [line.strip() for line in handle]
    with gzip.open(raw / "GSM5866081_features.tsv.gz", "rt") as handle:
        features = [line.rstrip().split("\t") for line in handle]

    feature_index = {
        index + 1: row[1]
        for index, row in enumerate(features)
        if row[1] in GENES
    }
    totals = [0] * len(cells)
    values = {gene: {} for gene in GENES}
    with gzip.open(raw / "GSM5866081_matrix.mtx.gz", "rt") as handle:
        line = next(handle)
        if not line.startswith("%%MatrixMarket matrix coordinate"):
            raise ValueError("Expected coordinate MTX")
        line = next(handle)
        while line.startswith("%"):
            line = next(handle)
        n_features, n_cells, n_nonzero = map(int, line.split())
        if (n_features, n_cells) != (len(features), len(cells)):
            raise ValueError("RNA matrix dimensions do not match features/barcodes")
        seen = 0
        for line in handle:
            feature, cell, value = map(int, line.split())
            cell -= 1
            totals[cell] += value
            gene = feature_index.get(feature)
            if gene is not None:
                values[gene][cell] = value
            seen += 1
        if seen != n_nonzero:
            raise ValueError("RNA matrix nonzero count is incomplete")

    groups = defaultdict(list)
    for index, cell_id in enumerate(cells):
        sample, barcode = cell_id.split("_", 1)
        record = audit.get(barcode)
        if record and record["n_matching_RNA_samples"] == "1":
            if sample != record["matching_RNA_sample_prefixes"]:
                raise ValueError("Unique barcode/sample mapping mismatch")
            donor = samples[sample]["donor"]
            groups[f"MG_{donor}"].append(index)
            groups["MG_all"].append(index)
        elif record:
            groups["ambiguous_excluded"].append(index)
        else:
            groups["non_MG_barcode_background"].append(index)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["group", "gene", "n_cells", "fraction_detected", "mean_log1p_CP10k"]
        )
        for group, indices in sorted(groups.items()):
            for gene in GENES:
                observed = values[gene]
                writer.writerow(
                    [
                        group,
                        gene,
                        len(indices),
                        sum(observed.get(i, 0) > 0 for i in indices) / len(indices),
                        sum(
                            math.log1p(observed.get(i, 0) * 10000 / totals[i])
                            if totals[i]
                            else 0
                            for i in indices
                        )
                        / len(indices),
                    ]
                )
    print(output)


if __name__ == "__main__":
    run()
