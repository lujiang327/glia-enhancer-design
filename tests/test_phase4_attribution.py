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


if __name__ == "__main__":
    unittest.main()
