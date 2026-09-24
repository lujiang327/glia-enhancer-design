#!/usr/bin/env python3
"""Combine cell-type fragment audits into a donor-by-cell-type QC table."""

import argparse
import csv
import json
from pathlib import Path


DONORS = ("LGS1", "LGS2", "LGS3", "LVG1")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-root", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args()

    summaries = []
    for path in sorted(args.audit_root.glob("*/audit_summary.json")):
        summaries.append(json.loads(path.read_text()))
    if len(summaries) != 13:
        raise ValueError("Expected 13 cell-type summaries, found {}".format(len(summaries)))

    rows = []
    cautions = []
    for summary in sorted(summaries, key=lambda x: x["analysis_group"]):
        for donor in DONORS:
            values = summary["donors"].get(donor, {})
            row = {
                "analysis_group": summary["analysis_group"],
                "published_label": summary["published_label"],
                "donor": donor,
                "unambiguous_cells": int(values.get("unambiguous_cells", 0)),
                "fragment_records": int(values.get("fragment_records", 0)),
                "canonical_fragment_records": int(values.get("canonical_fragment_records", 0)),
                "read_support_sum": int(values.get("read_support_sum", 0)),
                "canonical_read_support_sum": int(values.get("canonical_read_support_sum", 0)),
            }
            rows.append(row)
            if row["unambiguous_cells"] == 0:
                cautions.append("{} has no unambiguous cells for {}".format(donor, summary["analysis_group"]))
            elif row["unambiguous_cells"] < 50:
                cautions.append("{} has fewer than 50 cells for {}".format(donor, summary["analysis_group"]))

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    result = {
        "status": "CELL_TYPE_DONOR_AUDIT_READY_FOR_REVIEW",
        "cell_types": len(summaries),
        "donors": list(DONORS),
        "rows": rows,
        "cautions": cautions,
        "next_gate": (
            "Review donor coverage before creating donor-by-cell-type pseudobulks "
            "and calling reproducible peaks."
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
