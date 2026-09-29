import csv
import json
from pathlib import Path
import unittest


class Phase3DifferentialAccessibilityContract(unittest.TestCase):
    def test_primary_and_pairwise_contrasts_are_balanced(self):
        root = Path(__file__).resolve().parents[1]
        with (root / "config/phase3_cell_types.tsv").open(newline="") as handle:
            groups = [row["analysis_group"] for row in csv.DictReader(handle, delimiter="\t")]
        self.assertEqual(groups[0], "MG")
        self.assertEqual(len(groups), 13)
        off_targets = groups[1:]
        primary = {group: -1 / len(off_targets) for group in off_targets}
        primary["MG"] = 1
        self.assertAlmostEqual(sum(primary.values()), 0)
        self.assertEqual(len(set(primary[group] for group in off_targets)), 1)
        for group in off_targets:
            contrast = {value: 0 for value in groups}
            contrast["MG"] = 1
            contrast[group] = -1
            self.assertEqual(sum(contrast.values()), 0)

    def test_candidate_selection_remains_blocked(self):
        root = Path(__file__).resolve().parents[1]
        contract = json.loads((root / "config/phase3_differential_accessibility.json").read_text())
        self.assertEqual(
            contract["candidate_selection_status"],
            "BLOCKED_PENDING_DIFFERENTIAL_MODEL_QC",
        )


if __name__ == "__main__":
    unittest.main()
