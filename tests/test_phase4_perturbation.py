from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_phase4_perturbation_pilot import edit_screen_pass, load_pilot_candidates, nominate_edits, reverse_complement_one_hot, select_spaced


class Phase4Perturbation(unittest.TestCase):
    def test_full_scan_selects_all_eligible_parents(self):
        root = Path(__file__).resolve().parents[1]
        parents = load_pilot_candidates(
            root / "reports/phase4/attribution/attribution_qc.tsv.gz",
            root / "reports/phase4/parent_scoring/parent_predictions.tsv.gz", 947,
        )
        self.assertEqual(len(parents), 947)
        self.assertEqual(parents.peak_id.nunique(), 947)
        self.assertTrue(parents.attribution_qc_pass.all())

    def test_screen_requires_effect_and_both_profile_limits(self):
        effect = np.array([.1, -.1, .2, .2, .099])
        forward = np.array([.05, 0, .051, 0, 0])
        reverse = np.array([.05, 0, 0, .051, 0])
        self.assertEqual(edit_screen_pass(effect, forward, reverse, .1, .05).tolist(),
                         [True, False, False, False, False])

    def test_real_qc_joins_verified_parent_sequences(self):
        root = Path(__file__).resolve().parents[1]
        pilot = load_pilot_candidates(
            root / "reports/phase4/attribution/attribution_qc.tsv.gz",
            root / "reports/phase4/parent_scoring/parent_predictions.tsv.gz", 100,
        )
        self.assertEqual(len(pilot), 100)
        self.assertTrue(pilot.attribution_qc_pass.all())
        self.assertTrue((pilot.parent_sequence.str.len() == 500).all())
        self.assertEqual(pilot.candidate_rank.max(), 104)

    def test_select_spaced_prefers_score_then_position(self):
        records = [
            {"nomination_score": 3, "parent_position_0based": 10},
            {"nomination_score": 2, "parent_position_0based": 12},
            {"nomination_score": 1, "parent_position_0based": 20},
        ]
        selected = select_spaced(records, 2, 5)
        self.assertEqual([value["parent_position_0based"] for value in selected], [10, 20])

    def test_nomination_uses_only_parent_slice_and_alt_alleles(self):
        sequence = np.zeros((8, 4), dtype=np.int8)
        sequence[:, 0] = 1
        counts = np.zeros((8, 4), dtype=float)
        profile = np.zeros((8, 4), dtype=float)
        for position in range(2, 6):
            counts[position] = [0, position, -position, 0.5]
        edits = nominate_edits(sequence, counts, profile, 2, 6, 2, 2, 1)
        self.assertEqual(len(edits), 4)
        self.assertTrue(all(0 <= value["parent_position_0based"] < 4 for value in edits))
        self.assertTrue(all(value["ref"] != value["alt"] for value in edits))
        self.assertEqual({value["design_class"] for value in edits}, {"gain", "loss_control"})

    def test_reverse_complement_round_trip(self):
        values = np.arange(2 * 7 * 4).reshape(2, 7, 4)
        self.assertTrue(np.array_equal(reverse_complement_one_hot(reverse_complement_one_hot(values)), values))


if __name__ == "__main__":
    unittest.main()
