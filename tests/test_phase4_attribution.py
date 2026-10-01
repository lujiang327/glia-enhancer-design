from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from attribute_phase4_parent_sequences import (
    align_reverse_complement_scores,
    reverse_complement_one_hot,
    row_agreement,
)
from review_phase4_attribution import failure_reasons


class Phase4Attribution(unittest.TestCase):
    def test_reverse_complement_round_trip(self):
        values = np.arange(2 * 5 * 4).reshape(2, 5, 4)
        self.assertTrue(np.array_equal(reverse_complement_one_hot(reverse_complement_one_hot(values)), values))
        self.assertTrue(np.array_equal(align_reverse_complement_scores(align_reverse_complement_scores(values)), values))

    def test_identical_projected_maps_have_perfect_agreement(self):
        one_hot = np.zeros((1, 6, 4), dtype=float)
        one_hot[0, :, 0] = 1
        scores = np.zeros_like(one_hot)
        scores[0, :, 0] = [-2, -1, 0, 1, 2, 3]
        projected = scores * one_hot
        pearson, cosine, sign = row_agreement(projected, projected, 0, 6)
        self.assertAlmostEqual(pearson[0], 1)
        self.assertAlmostEqual(cosine[0], 1)
        self.assertAlmostEqual(sign[0], 1)

    def test_review_gate_records_specific_failures(self):
        row = {
            "count_attribution_parent_pearson": 0.8,
            "count_attribution_parent_cosine": 0.7,
            "count_attribution_parent_top10pct_sign_concordance": 0.9,
            "profile_attribution_parent_pearson": 0.4,
            "profile_attribution_parent_cosine": 0.3,
            "profile_attribution_parent_top10pct_sign_concordance": 0.7,
            "model_input_invalid_bases": 2,
        }
        self.assertEqual(
            failure_reasons(row),
            "profile_pearson_lt_0.5;profile_cosine_lt_0.5;"
            "profile_top10pct_sign_lt_0.8;ambiguous_model_input_bases",
        )


if __name__ == "__main__":
    unittest.main()
