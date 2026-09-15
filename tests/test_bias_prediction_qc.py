import json
from pathlib import Path
import tempfile
import unittest

from scripts.summarize_bias_prediction_qc import run


def metrics(nonpeak_pearson, peak_pearson):
    result = {'counts_metrics': {}, 'profile_metrics': {}}
    for group, pearson in [('nonpeaks', nonpeak_pearson), ('peaks', peak_pearson)]:
        result['counts_metrics'][group] = {'pearsonr': pearson, 'spearmanr': 0.2, 'mse': 0.4}
        result['profile_metrics'][group] = {'median_jsd': 0.3, 'median_norm_jsd': 0.5}
    return result


class BiasPredictionQCTest(unittest.TestCase):
    def summarize(self, payload):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'metrics.json'
            output = Path(directory) / 'summary.json'
            source.write_text(json.dumps(payload))
            run(source, output)
            return json.loads(output.read_text())

    def test_pass_requires_positive_nonpeak_and_acceptable_peak_correlation(self):
        result = self.summarize(metrics(0.2, -0.1))
        self.assertEqual(result['status'], 'PASS_NUMERICAL_QC_MOTIF_QC_PENDING')

    def test_peak_caution_range_is_preserved(self):
        result = self.summarize(metrics(0.2, -0.4))
        self.assertEqual(result['status'], 'CAUTION_NUMERICAL_QC_MOTIF_QC_REQUIRED')

    def test_nonpositive_nonpeak_correlation_fails(self):
        result = self.summarize(metrics(0.0, -0.1))
        self.assertEqual(result['status'], 'FAIL_NUMERICAL_BIAS_QC')


if __name__ == '__main__':
    unittest.main()
