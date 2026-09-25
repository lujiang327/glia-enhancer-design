#!/usr/bin/env python3
"""Combine 13 Phase 3 pseudobulk summaries into a reviewable QC table."""

import argparse
import csv
import json
import statistics
from pathlib import Path


DONORS = ("LGS1", "LGS2", "LGS3", "LVG1")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pseudobulk-root", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args()

    paths = sorted(args.pseudobulk_root.glob("*/pseudobulk_summary.json"))
    if len(paths) != 13:
        raise ValueError("Expected 13 pseudobulk summaries, found {}".format(len(paths)))
    summaries = [json.loads(path.read_text()) for path in paths]
    groups = [summary["analysis_group"] for summary in summaries]
    if len(set(groups)) != 13:
        raise ValueError("Duplicate analysis groups in pseudobulk summaries")

    rows = []
    cautions = []
    for summary in sorted(summaries, key=lambda value: value["analysis_group"]):
        depths = []
        for donor in DONORS:
            values = summary["donors"][donor]
            row = {
                "analysis_group": summary["analysis_group"],
                "published_label": summary["published_label"],
                "donor": donor,
                "unambiguous_cells": int(values["unambiguous_cells"]),
                "retained_fragment_records": int(values["retained_fragment_records"]),
                "retained_read_support": int(values["retained_read_support"]),
                "blacklist_overlap_records": int(values.get("blacklist_overlap", 0)),
                "unique_to_support_ratio": float(values["unique_to_support_ratio"]),
                "output_bytes": int(values["output_bytes"]),
                "output_sha256": values["output_sha256"],
            }
            rows.append(row)
            depths.append(row["retained_fragment_records"])
            if row["retained_fragment_records"] < 250000:
                cautions.append(
                    "{} {} has fewer than 250,000 retained fragments".format(
                        summary["analysis_group"], donor
                    )
                )
        mean_depth = statistics.mean(depths)
        cv = statistics.pstdev(depths) / mean_depth if mean_depth else None
        if cv is not None and cv > 0.75:
            cautions.append(
                "{} donor depth CV exceeds 0.75 ({:.3f})".format(
                    summary["analysis_group"], cv
                )
            )

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    result = {
        "status": "PSEUDOBULK_QC_READY_FOR_REVIEW",
        "cell_types": len(summaries),
        "pseudobulks": len(rows),
        "donors": list(DONORS),
        "total_retained_fragment_records": sum(
            row["retained_fragment_records"] for row in rows
        ),
        "total_retained_read_support": sum(
            row["retained_read_support"] for row in rows
        ),
        "rows": rows,
        "cautions": cautions,
        "next_gate": (
            "Review all 52 donor pseudobulks before pooled and donor-level peak "
            "calling. Differential accessibility remains blocked."
        ),
    }
    args.output_json.write_text(json.dumps(result, indent=2) + "\n")
    with args.output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
