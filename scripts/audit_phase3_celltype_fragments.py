#!/usr/bin/env python3
"""Audit one published cell-type fragment file and assign unambiguous donors."""

import argparse
import collections
import csv
import gzip
import json
from pathlib import Path


CANONICAL = {"chr{}".format(i) for i in range(1, 23)} | {"chrX"}


def load_rna_assignments(barcodes_path):
    assignments = collections.defaultdict(set)
    with gzip.open(str(barcodes_path), "rt") as handle:
        for line in handle:
            cell_id = line.strip()
            if not cell_id:
                continue
            prefix, barcode = cell_id.split("_", 1)
            assignments[barcode].add(prefix)
    return assignments


def load_prefix_donors(metadata_path):
    result = {}
    with metadata_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            prefix, donor = row["rna_prefix"], row["donor"]
            if prefix in result and result[prefix] != donor:
                raise ValueError("RNA prefix maps to multiple donors: {}".format(prefix))
            result[prefix] = donor
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fragments", type=Path, required=True)
    parser.add_argument("--rna-barcodes", type=Path, required=True)
    parser.add_argument("--donor-metadata", type=Path, required=True)
    parser.add_argument("--file-label", required=True)
    parser.add_argument("--published-label", required=True)
    parser.add_argument("--analysis-group", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if not args.fragments.is_file():
        raise FileNotFoundError(args.fragments)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    assignments = load_rna_assignments(args.rna_barcodes)
    prefix_donors = load_prefix_donors(args.donor_metadata)

    record_counts = collections.Counter()
    canonical_counts = collections.Counter()
    read_support = collections.Counter()
    canonical_read_support = collections.Counter()
    total_records = 0
    malformed = 0
    with gzip.open(str(args.fragments), "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                malformed += 1
                continue
            chrom, start, end, barcode, support = fields
            start, end, support = int(start), int(end), int(support)
            if start < 0 or end <= start or support < 1:
                malformed += 1
                continue
            total_records += 1
            record_counts[barcode] += 1
            read_support[barcode] += support
            if chrom in CANONICAL:
                canonical_counts[barcode] += 1
                canonical_read_support[barcode] += support
            if total_records % 10000000 == 0:
                print("{}: {:,} records".format(args.analysis_group, total_records), flush=True)
    if malformed:
        raise ValueError("{} malformed/invalid records".format(malformed))

    mapping_barcode_counts = collections.Counter()
    mapping_record_counts = collections.Counter()
    mapping_read_support = collections.Counter()
    donor_cells = collections.Counter()
    donor_records = collections.Counter()
    donor_canonical_records = collections.Counter()
    donor_read_support = collections.Counter()
    donor_canonical_read_support = collections.Counter()
    assignment_path = args.output_dir / "cell_assignments.tsv.gz"
    with gzip.open(str(assignment_path), "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "analysis_group", "published_label", "bare_barcode", "rna_prefix",
            "donor", "fragment_records", "canonical_fragment_records", "read_support_sum"
        ])
        for barcode in sorted(record_counts):
            prefixes = sorted(assignments.get(barcode, set()))
            mapping_class = str(len(prefixes))
            mapping_barcode_counts[mapping_class] += 1
            mapping_record_counts[mapping_class] += record_counts[barcode]
            mapping_read_support[mapping_class] += read_support[barcode]
            if len(prefixes) != 1:
                continue
            prefix = prefixes[0]
            donor = prefix_donors.get(prefix)
            if donor is None:
                raise ValueError("No donor for RNA prefix {}".format(prefix))
            donor_cells[donor] += 1
            donor_records[donor] += record_counts[barcode]
            donor_canonical_records[donor] += canonical_counts[barcode]
            donor_read_support[donor] += read_support[barcode]
            donor_canonical_read_support[donor] += canonical_read_support[barcode]
            writer.writerow([
                args.analysis_group, args.published_label, barcode, prefix, donor,
                record_counts[barcode], canonical_counts[barcode], read_support[barcode]
            ])

    donors = sorted(set(prefix_donors.values()))
    summary = {
        "file_label": args.file_label,
        "published_label": args.published_label,
        "analysis_group": args.analysis_group,
        "fragment_file": str(args.fragments),
        "fragment_records": total_records,
        "read_support_sum": sum(read_support.values()),
        "distinct_bare_barcodes": len(record_counts),
        "mapping_barcode_counts": dict(mapping_barcode_counts),
        "mapping_record_counts": dict(mapping_record_counts),
        "mapping_read_support": dict(mapping_read_support),
        "donors": {
            donor: {
                "unambiguous_cells": donor_cells[donor],
                "fragment_records": donor_records[donor],
                "canonical_fragment_records": donor_canonical_records[donor],
                "read_support_sum": donor_read_support[donor],
                "canonical_read_support_sum": donor_canonical_read_support[donor],
            }
            for donor in donors
        },
        "policy": (
            "Published cell-type file labels are retained. A bare barcode is assigned "
            "only when it matches exactly one prefixed RNA barcode; ambiguous barcodes "
            "are excluded from donor-level pseudobulks."
        ),
    }
    (args.output_dir / "audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
