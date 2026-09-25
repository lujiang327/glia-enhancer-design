#!/usr/bin/env python3
"""Build a fixed-width, donor-reproducible peak universe across cell types."""

import argparse
import bisect
import collections
import csv
import gzip
import json
from pathlib import Path


DONORS = ("LGS1", "LGS2", "LGS3", "LVG1")


def load_sizes(path):
    sizes = collections.OrderedDict()
    with path.open() as handle:
        for line in handle:
            chrom, size = line.split()[:2]
            sizes[chrom] = int(size)
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
    return merge_intervals(intervals)


def merge_intervals(intervals):
    result = {}
    for chrom, values in intervals.items():
        merged = []
        for start, end in sorted(values):
            if current_overlap(merged, start):
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        result[chrom] = (
            [value[0] for value in merged],
            [value[1] for value in merged],
        )
    return result


def current_overlap(merged, start):
    return bool(merged and start <= merged[-1][1])


def overlaps(index, chrom, start, end):
    starts, ends = index.get(chrom, ((), ()))
    position = bisect.bisect_left(starts, end) - 1
    return position >= 0 and ends[position] > start


def load_peak_rows(path):
    rows = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                raise ValueError("Expected narrowPeak columns at {}:{}".format(path, line_number))
            rows.append(fields)
    return rows


def donor_peak_index(path, sizes, blacklist):
    intervals = collections.defaultdict(list)
    for fields in load_peak_rows(path):
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        if chrom not in sizes or start < 0 or end <= start or end > sizes[chrom]:
            continue
        if overlaps(blacklist, chrom, start, end):
            continue
        intervals[chrom].append((start, end))
    return merge_intervals(intervals)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peak-root", type=Path, required=True)
    parser.add_argument("--cell-types", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--blacklist", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--width", type=int, default=500)
    parser.add_argument("--min-donors", type=int, default=2)
    args = parser.parse_args()

    if args.width < 2 or args.width % 2:
        raise ValueError("Peak width must be a positive even integer")
    if args.min_donors < 1 or args.min_donors > len(DONORS):
        raise ValueError("min-donors must be between 1 and 4")
    if args.output_dir.exists():
        raise FileExistsError("Refusing to overwrite: {}".format(args.output_dir))
    args.output_dir.mkdir(parents=True)

    sizes = load_sizes(args.chrom_sizes)
    chrom_rank = {chrom: index for index, chrom in enumerate(sizes)}
    blacklist = load_blacklist(args.blacklist)
    with args.cell_types.open(newline="") as handle:
        groups = [row["analysis_group"] for row in csv.DictReader(handle, delimiter="\t")]
    if len(groups) != 13 or len(set(groups)) != 13:
        raise ValueError("Expected 13 unique analysis groups")

    candidates = []
    group_qc = {}
    half_width = args.width // 2
    for group in groups:
        pooled_path = args.peak_root / group / "pooled" / "{}_pooled_peaks.narrowPeak".format(group)
        if not pooled_path.is_file():
            raise FileNotFoundError(pooled_path)
        donor_indexes = {}
        for donor in DONORS:
            path = args.peak_root / group / "donor_{}".format(donor) / "{}_{}_peaks.narrowPeak".format(group, donor)
            if not path.is_file():
                raise FileNotFoundError(path)
            donor_indexes[donor] = donor_peak_index(path, sizes, blacklist)

        counters = collections.Counter()
        for fields in load_peak_rows(pooled_path):
            counters["pooled_raw"] += 1
            chrom = fields[0]
            peak_start, summit_offset = int(fields[1]), int(fields[9])
            if chrom not in sizes or summit_offset < 0:
                counters["invalid_or_noncanonical"] += 1
                continue
            summit = peak_start + summit_offset
            start, end = summit - half_width, summit + half_width
            if start < 0 or end > sizes[chrom]:
                counters["boundary_excluded"] += 1
                continue
            if overlaps(blacklist, chrom, start, end):
                counters["blacklist_excluded"] += 1
                continue
            supporting = [
                donor
                for donor in DONORS
                if overlaps(donor_indexes[donor], chrom, start, end)
            ]
            if len(supporting) < args.min_donors:
                counters["insufficient_donor_support"] += 1
                continue
            signal = float(fields[6])
            pvalue = float(fields[7])
            qvalue = float(fields[8])
            candidates.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "summit": summit,
                    "source_group": group,
                    "donors": supporting,
                    "donor_support": len(supporting),
                    "signal": signal,
                    "pvalue": pvalue,
                    "qvalue": qvalue,
                }
            )
            counters["reproducible_candidates"] += 1
        group_qc[group] = dict(counters)

    # Resolve overlaps using donor reproducibility first, then MACS3 strength.
    accepted_index = {chrom: ([], []) for chrom in sizes}
    accepted = []
    ordered = sorted(
        candidates,
        key=lambda row: (
            -row["donor_support"],
            -row["qvalue"],
            -row["signal"],
            row["source_group"],
            chrom_rank[row["chrom"]],
            row["start"],
        ),
    )
    for candidate in ordered:
        starts, ends = accepted_index[candidate["chrom"]]
        position = bisect.bisect_left(starts, candidate["start"])
        overlaps_previous = position > 0 and ends[position - 1] > candidate["start"]
        overlaps_next = position < len(starts) and starts[position] < candidate["end"]
        if overlaps_previous or overlaps_next:
            continue
        starts.insert(position, candidate["start"])
        ends.insert(position, candidate["end"])
        accepted.append(candidate)

    accepted.sort(key=lambda row: (chrom_rank[row["chrom"]], row["start"], row["end"]))
    bed_path = args.output_dir / "retina_consensus_peaks.bed"
    narrow_path = args.output_dir / "retina_consensus_peaks.narrowPeak"
    metadata_path = args.output_dir / "retina_consensus_peak_metadata.tsv"
    with bed_path.open("w") as bed, narrow_path.open("w") as narrow, metadata_path.open("w", newline="") as meta:
        fieldnames = [
            "peak_id", "chrom", "start", "end", "source_group", "donor_support",
            "supporting_donors", "macs_signal", "macs_pvalue", "macs_qvalue",
        ]
        writer = csv.DictWriter(meta, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for index, row in enumerate(accepted, 1):
            peak_id = "retina_peak_{:06d}".format(index)
            score = min(1000, max(0, int(round(row["qvalue"] * 10))))
            bed.write("{}\t{}\t{}\t{}\n".format(row["chrom"], row["start"], row["end"], peak_id))
            narrow.write(
                "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\t{}\n".format(
                    row["chrom"], row["start"], row["end"], peak_id, score,
                    row["signal"], row["pvalue"], row["qvalue"], half_width,
                )
            )
            writer.writerow(
                {
                    "peak_id": peak_id,
                    "chrom": row["chrom"],
                    "start": row["start"],
                    "end": row["end"],
                    "source_group": row["source_group"],
                    "donor_support": row["donor_support"],
                    "supporting_donors": ",".join(row["donors"]),
                    "macs_signal": row["signal"],
                    "macs_pvalue": row["pvalue"],
                    "macs_qvalue": row["qvalue"],
                }
            )

    selected_by_group = collections.Counter(row["source_group"] for row in accepted)
    result = {
        "status": "CONSENSUS_PEAK_UNIVERSE_READY_FOR_COUNTING",
        "width": args.width,
        "minimum_donor_support": args.min_donors,
        "cell_types": groups,
        "group_qc": group_qc,
        "reproducible_candidates_before_overlap_resolution": len(candidates),
        "consensus_peaks": len(accepted),
        "selected_source_counts": dict(sorted(selected_by_group.items())),
        "overlap_policy": (
            "Fixed-width peaks are greedily selected by donor support, MACS3 q-value, "
            "signal, source group, chromosome, and coordinate; final peaks do not overlap."
        ),
        "outputs": {
            "bed": str(bed_path),
            "narrowPeak": str(narrow_path),
            "metadata": str(metadata_path),
        },
    }
    (args.output_dir / "consensus_peak_summary.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
