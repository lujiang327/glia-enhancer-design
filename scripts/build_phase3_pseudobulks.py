#!/usr/bin/env python3
"""Build four donor pseudobulks from one published cell-type fragment pool.

Only barcodes accepted by the preceding unambiguous RNA-prefix audit are used.
Each retained fragment row is written once as three-column BEDPE. Column-five
read support is summarized for QC but is not used to duplicate molecules.
"""

import argparse
import bisect
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path


DONORS = ("LGS1", "LGS2", "LGS3", "LVG1")


def load_sizes(path):
    sizes = collections.OrderedDict()
    with path.open() as handle:
        for line in handle:
            chrom, size = line.split()[:2]
            sizes[chrom] = int(size)
    expected = {"chr{}".format(i) for i in range(1, 23)} | {"chrX"}
    if set(sizes) != expected:
        raise ValueError("Chromosome sizes must contain exactly chr1-22 and chrX")
    return sizes


def load_blacklist(path):
    intervals = collections.defaultdict(list)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(str(path), "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            chrom, start, end = line.split()[:3]
            intervals[chrom].append((int(start), int(end)))
    merged = {}
    for chrom, values in intervals.items():
        current = []
        for start, end in sorted(values):
            if current and start <= current[-1][1]:
                current[-1] = (current[-1][0], max(current[-1][1], end))
            else:
                current.append((start, end))
        merged[chrom] = (
            [value[0] for value in current],
            [value[1] for value in current],
        )
    return merged


def overlaps(intervals, chrom, start, end):
    starts, ends = intervals.get(chrom, ((), ()))
    index = bisect.bisect_left(starts, end) - 1
    return index >= 0 and ends[index] > start


def load_assignments(path, expected_group):
    assignments = {}
    expected = {donor: collections.Counter() for donor in DONORS}
    with gzip.open(str(path), "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["analysis_group"] != expected_group:
                raise ValueError("Unexpected analysis group in assignment table")
            barcode, donor = row["bare_barcode"], row["donor"]
            if donor not in expected:
                raise ValueError("Unexpected donor: {}".format(donor))
            if barcode in assignments:
                raise ValueError("Duplicate accepted barcode: {}".format(barcode))
            assignments[barcode] = donor
            expected[donor]["cells"] += 1
            expected[donor]["fragment_records"] += int(row["fragment_records"])
            expected[donor]["canonical_fragment_records"] += int(
                row["canonical_fragment_records"]
            )
            expected[donor]["read_support_sum"] += int(row["read_support_sum"])
    if not assignments:
        raise ValueError("Assignment table is empty")
    return assignments, expected


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(2 ** 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fragments", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--blacklist", type=Path, required=True)
    parser.add_argument("--analysis-group", required=True)
    parser.add_argument("--published-label", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    for path in (args.fragments, args.assignments, args.chrom_sizes, args.blacklist):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.output_dir.exists():
        raise FileExistsError("Refusing to overwrite: {}".format(args.output_dir))
    args.output_dir.mkdir(parents=True)

    sizes = load_sizes(args.chrom_sizes)
    blacklist = load_blacklist(args.blacklist)
    assignments, expected = load_assignments(args.assignments, args.analysis_group)
    flow = collections.Counter()
    donor_counts = {donor: collections.Counter() for donor in DONORS}
    retained_cells = {donor: set() for donor in DONORS}
    part_paths = {
        donor: args.output_dir / "{}.fragments.bedpe.gz.part".format(donor)
        for donor in DONORS
    }
    final_paths = {
        donor: args.output_dir / "{}.fragments.bedpe.gz".format(donor)
        for donor in DONORS
    }
    handles = {
        donor: gzip.open(str(part_paths[donor]), "wt", compresslevel=1)
        for donor in DONORS
    }
    try:
        with gzip.open(str(args.fragments), "rt") as source:
            for line_number, line in enumerate(source, 1):
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 5:
                    raise ValueError("Expected five columns at line {}".format(line_number))
                chrom, start_text, end_text, barcode, support_text = fields
                start, end, support = int(start_text), int(end_text), int(support_text)
                flow["input_records"] += 1
                flow["input_read_support"] += support
                donor = assignments.get(barcode)
                if donor is None:
                    flow["barcode_not_unambiguous"] += 1
                    flow["barcode_not_unambiguous_support"] += support
                    continue
                donor_counts[donor]["assigned_records"] += 1
                donor_counts[donor]["assigned_read_support"] += support
                if chrom not in sizes:
                    flow["noncanonical"] += 1
                    continue
                if start < 0 or end <= start or end > sizes[chrom] or support < 1:
                    flow["invalid"] += 1
                    continue
                donor_counts[donor]["canonical_records_before_blacklist"] += 1
                donor_counts[donor]["canonical_support_before_blacklist"] += support
                if overlaps(blacklist, chrom, start, end):
                    flow["blacklist_overlap"] += 1
                    donor_counts[donor]["blacklist_overlap"] += 1
                    donor_counts[donor]["blacklist_support"] += support
                    continue
                handles[donor].write("{}\t{}\t{}\n".format(chrom, start, end))
                retained_cells[donor].add(barcode)
                donor_counts[donor]["retained_fragment_records"] += 1
                donor_counts[donor]["retained_read_support"] += support
                flow["retained_fragment_records"] += 1
                flow["retained_read_support"] += support
                if flow["input_records"] % 10000000 == 0:
                    print(
                        "{}: {:,} input; {:,} retained".format(
                            args.analysis_group,
                            flow["input_records"],
                            flow["retained_fragment_records"],
                        ),
                        flush=True,
                    )
    finally:
        for handle in handles.values():
            handle.close()

    for donor in DONORS:
        observed = donor_counts[donor]
        if observed["assigned_records"] != expected[donor]["fragment_records"]:
            raise ValueError("Assigned-record mismatch for {}".format(donor))
        if (
            observed["canonical_records_before_blacklist"]
            != expected[donor]["canonical_fragment_records"]
        ):
            raise ValueError("Canonical-record mismatch for {}".format(donor))
        if observed["assigned_read_support"] != expected[donor]["read_support_sum"]:
            raise ValueError("Read-support mismatch for {}".format(donor))
        if len(retained_cells[donor]) != expected[donor]["cells"]:
            raise ValueError("Retained-cell mismatch for {}".format(donor))
        part_paths[donor].replace(final_paths[donor])

    donors = {}
    for donor in DONORS:
        values = dict(donor_counts[donor])
        values["unambiguous_cells"] = len(retained_cells[donor])
        values["output"] = str(final_paths[donor])
        values["output_bytes"] = final_paths[donor].stat().st_size
        values["output_sha256"] = sha256(final_paths[donor])
        values["unique_to_support_ratio"] = (
            values["retained_fragment_records"] / values["retained_read_support"]
            if values["retained_read_support"]
            else None
        )
        donors[donor] = values
    result = {
        "status": "PSEUDOBULKS_READY_FOR_QC_REVIEW",
        "analysis_group": args.analysis_group,
        "published_label": args.published_label,
        "fragment_file": str(args.fragments),
        "assignment_file": str(args.assignments),
        "coordinate_policy": "chr1-22/X, valid hg38 coordinates, full-fragment blacklist exclusion",
        "counting_policy": (
            "One retained fragment-file row is one molecule. Fifth-column read "
            "support is summarized but does not duplicate rows."
        ),
        "flow": dict(flow),
        "donors": donors,
    }
    (args.output_dir / "pseudobulk_summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
