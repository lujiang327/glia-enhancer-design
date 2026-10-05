from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from annotate_phase5_motif_attribution import project, site_features, edit_features


class MotifAttribution(unittest.TestCase):
    def test_projection_uses_observed_base_not_largest_channel(self):
        scores = np.array([[1,99,0,0], [99,2,0,0], [0,0,-3,99], [99,0,0,4]])
        observed = np.eye(4)
        self.assertEqual(project(scores, observed).tolist(), [1,2,-3,4])
        observed[0] = 0
        with self.assertRaises(ValueError):
            project(scores, observed)

    def test_site_half_open_and_opposite_orientation_signs(self):
        f = np.array([100.,2.,3.,-100.])
        r = np.array([-100.,-1.,-2.,100.])
        c = (f+r)/2
        result = site_features(f,r,c,1,3)
        self.assertEqual(result['site_sum_forward'],5)
        self.assertEqual(result['site_sum_reverse_complement_aligned'],-3)
        self.assertFalse(result['site_positive_both_orientations'])
        self.assertFalse(result['site_negative_both_orientations'])
        self.assertEqual(result['site_absolute_attribution_fraction_of_parent'],1)
        with self.assertRaises(ValueError):
            site_features(f,r,c,1,5)

    def test_alt_minus_ref_channel_order_and_ties(self):
        h = np.array([[2.,-1.,8.,3.], [0.,0.,0.,0.]])
        projected = np.array([2.,0.])
        result = edit_features(h,projected,0,'A','G')
        self.assertEqual(result['alt_minus_ref_hypothetical_delta'],6)
        self.assertEqual(result['edited_parent_base_attribution'],2)
        self.assertEqual(result['edited_base_absolute_attribution_percentile'],1)
        self.assertEqual(edit_features(h,projected,0,'G','A')['alt_minus_ref_hypothetical_delta'],-6)
        with self.assertRaises(ValueError):
            edit_features(h,projected,0,'A','A')

if __name__ == '__main__':
    unittest.main()
