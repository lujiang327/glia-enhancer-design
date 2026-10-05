from pathlib import Path
import sys
import unittest
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from annotate_phase5_motif_attribution import project, site_features, edit_features
from review_phase5_motif_attribution import contextual_flags


class MotifAttribution(unittest.TestCase):
    def test_enzyme_and_low_attribution_matches_are_separate_context(self):
        sites = pd.DataFrame(dict(motif_id=['TN5_1','DNASE_2','NFI','SOX9','AP1'],
            counts_site_positive_both_orientations=[True,True,True,True,False],
            counts_site_negative_both_orientations=[False,False,False,False,True],
            counts_site_top10pct_parent_bases=[2,2,2,0,1]))
        result = contextual_flags(sites)
        self.assertEqual(result.enzyme_motif_reference_match.tolist(), [True,True,False,False,False])
        self.assertEqual(result.nonenzyme_positive_parent_count_context.tolist(), [False,False,True,False,False])
        self.assertEqual(result.nonenzyme_negative_parent_count_context.tolist(), [False,False,False,False,True])

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
