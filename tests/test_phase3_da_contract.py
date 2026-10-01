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

    def test_candidate_ranking_gate_matches_review(self):
        root = Path(__file__).resolve().parents[1]
        contract = json.loads((root / "config/phase3_differential_accessibility.json").read_text())
        self.assertEqual(
            contract["candidate_selection_status"],
            "GO_FOR_CANDIDATE_RANKING_WITH_CAVEATS",
        )
        self.assertEqual(
            contract["qc_review"]["decision"],
            contract["candidate_selection_status"],
        )
        self.assertTrue(contract["qc_review"]["primary_eligibility"]["positive_in_all_four_donors"])


if __name__ == "__main__":
    unittest.main()
