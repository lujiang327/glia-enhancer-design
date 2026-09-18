import json
from pathlib import Path
import tempfile
import unittest

from scripts.compare_chrombpnet_prediction_qc import run


def metrics(pearson, spearman, mse, jsd, normalized_jsd):
    result = {"counts_metrics": {}, "profile_metrics": {}}
    for group in ("peaks", "nonpeaks", "peaks_and_nonpeaks"):
        result["counts_metrics"][group] = {
            "pearsonr": pearson,
            "spearmanr": spearman,
            "mse": mse,
        }
        result["profile_metrics"][group] = {
            "median_jsd": jsd,
            "median_norm_jsd": normalized_jsd,
        }
    return result


class ChromBPNetPredictionQCTest(unittest.TestCase):
    def summarize(self, model, bias):
        with tempfile.TemporaryDirectory() as directory:
            model_path = Path(directory) / "model.json"
            bias_path = Path(directory) / "bias.json"
            output = Path(directory) / "summary.json"
            model_path.write_text(json.dumps(model))
            bias_path.write_text(json.dumps(bias))
            run(model_path, bias_path, output)
            return json.loads(output.read_text())

    def test_positive_deltas_always_favor_full_model(self):
        model = metrics(0.8, 0.7, 0.5, 0.3, 0.6)
        bias = metrics(0.3, 0.2, 2.0, 0.6, 0.2)
        result = self.summarize(model, bias)
        peak = result["improvement_over_scaled_bias"]["peaks"]
        self.assertAlmostEqual(peak["counts_pearsonr_gain"], 0.5)
        self.assertAlmostEqual(peak["counts_mse_reduction"], 1.5)
        self.assertAlmostEqual(peak["median_jsd_reduction"], 0.3)
        self.assertEqual(result["cautions"], [])

    def test_peak_regressions_are_reported_for_review(self):
        model = metrics(0.2, 0.1, 2.0, 0.7, 0.1)
        bias = metrics(0.3, 0.2, 1.0, 0.6, 0.2)
        result = self.summarize(model, bias)
        self.assertEqual(len(result["cautions"]), 5)


if __name__ == "__main__":
    unittest.main()
