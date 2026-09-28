#!/usr/bin/env python3
"""Count one donor-by-cell-type pseudobulk on the Phase 3 peak universe.

Each retained fragment-file row contributes one count to every consensus peak
that it overlaps.  Consensus peaks are nonoverlapping, and fragments that span
more than one peak are counted explicitly in the QC summary.
"""

import argparse
import bisect
import collections
import gzip
import hashlib
import json
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(2 ** 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_peaks(path):
    peak_rows = collections.defaultdict(list)
    ordered = []
    previous = None
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise ValueError("Expected four BED columns at {}:{}".format(path, line_number))
            chrom, start_text, end_text, peak_id = fields
            start, end = int(start_text), int(end_text)
            if end <= start:
                raise ValueError("Invalid peak interval at {}:{}".format(path, line_number))
            if previous is not None and chrom == previous[0] and start < previous[2]:
                raise ValueError("Peak BED is overlapping or unsorted within chromosome")
            if peak_id != "retina_peak_{:06d}".format(line_number):
                raise ValueError("Unexpected peak ID at line {}: {}".format(line_number, peak_id))
            peak_index = len(ordered)
            peak_rows[chrom].append((start, end, peak_index))
            ordered.append((chrom, start, end, peak_id))
            previous = (chrom, start, end)
    if not ordered:
        raise ValueError("Peak BED is empty")
    peaks = {}
    for chrom, rows in peak_rows.items():
        peaks[chrom] = (
            [row[0] for row in rows],
            [row[1] for row in rows],
            [row[2] for row in rows],
        )
    return peaks, ordered


def count_fragments(fragment_path, peaks, peak_count):
    counts = [0] * peak_count
    flow = collections.Counter()
    with gzip.open(str(fragment_path), "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 3:
                raise ValueError(
                    "Expected three BEDPE columns at {}:{}".format(fragment_path, line_number)
                )
            chrom, start_text, end_text = fields
            start, end = int(start_text), int(end_text)
            if start < 0 or end <= start:
                raise ValueError("Invalid fragment at {}:{}".format(fragment_path, line_number))
            flow["fragment_records"] += 1
            starts, ends, indexes = peaks.get(chrom, ((), (), ()))
            # The released cell-type fragment pools are not coordinate sorted.
            # Find the first peak whose end is after this fragment's start;
            # bisect runs in C and keeps counting independent of input order.
            position = bisect.bisect_right(ends, start)
            overlaps = 0
            while position < len(starts) and starts[position] < end:
                counts[indexes[position]] += 1
                overlaps += 1
                position += 1
            if overlaps:
                flow["fragments_in_any_peak"] += 1
                flow["peak_overlap_assignments"] += overlaps
                if overlaps > 1:
                    flow["fragments_overlapping_multiple_peaks"] += 1
                    flow["extra_peak_assignments"] += overlaps - 1
            if line_number % 10000000 == 0:
                print(
                    "{:,} fragments; {:,} in peaks".format(
                        flow["fragment_records"], flow["fragments_in_any_peak"]
                    ),
                    flush=True,
                )
    return counts, flow


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fragments", type=Path, required=True)
    parser.add_argument("--peaks", type=Path, required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--analysis-group", required=True)
    parser.add_argument("--donor", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    for path in (args.fragments, args.peaks):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.output_dir.exists():
        raise FileExistsError("Refusing to overwrite: {}".format(args.output_dir))
    args.output_dir.mkdir(parents=True)

    peaks, ordered = load_peaks(args.peaks)
    counts, flow = count_fragments(args.fragments, peaks, len(ordered))
    count_path = args.output_dir / "peak_counts.tsv.gz"
    with gzip.open(str(count_path) + ".part", "wt", compresslevel=1) as handle:
        handle.write("peak_id\tcount\n")
        for row, count in zip(ordered, counts):
            handle.write("{}\t{}\n".format(row[3], count))
    Path(str(count_path) + ".part").replace(count_path)

    fragments = flow["fragment_records"]
    in_peaks = flow["fragments_in_any_peak"]
    result = {
        "status": "PEAK_COUNTS_READY_FOR_MATRIX_QC",
        "sample_id": args.sample_id,
        "analysis_group": args.analysis_group,
        "donor": args.donor,
        "fragment_file": str(args.fragments),
        "peak_file": str(args.peaks),
        "peak_file_sha256": sha256(args.peaks),
        "counting_policy": (
            "Each retained fragment-file row contributes one count to every "
            "nonoverlapping consensus peak it intersects."
        ),
        "fragment_records": fragments,
        "fragments_in_any_peak": in_peaks,
        "frip": in_peaks / fragments if fragments else None,
        "peak_overlap_assignments": flow["peak_overlap_assignments"],
        "fragments_overlapping_multiple_peaks": flow[
            "fragments_overlapping_multiple_peaks"
        ],
        "extra_peak_assignments": flow["extra_peak_assignments"],
        "consensus_peaks": len(ordered),
        "detected_peaks": sum(value > 0 for value in counts),
        "zero_count_peaks": sum(value == 0 for value in counts),
        "maximum_peak_count": max(counts),
        "count_file": str(count_path),
        "count_file_sha256": sha256(count_path),
    }
    summary_path = args.output_dir / "count_summary.json"
    summary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    (args.output_dir / "output_checksums.sha256").write_text(
        "{}  {}\n{}  {}\n".format(
            sha256(count_path), count_path.name, sha256(summary_path), summary_path.name
        )
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
