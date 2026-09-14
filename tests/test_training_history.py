import csv
import json
from pathlib import Path
import tempfile
import unittest

from scripts.summarize_keras_history import run


class TrainingHistoryTest(unittest.TestCase):
    def test_summarizes_best_validation_epoch(self):
        with tempfile.TemporaryDirectory() as directory:
            history = Path(directory) / 'history.csv'
            output = Path(directory) / 'summary.json'
            with history.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=['epoch', 'loss', 'val_loss'])
                writer.writeheader()
                writer.writerows([
                    {'epoch': 0, 'loss': 3.0, 'val_loss': 4.0},
                    {'epoch': 1, 'loss': 2.0, 'val_loss': 2.5},
                    {'epoch': 2, 'loss': 1.5, 'val_loss': 3.0},
                ])
            run(history, output)
            result = json.loads(output.read_text())
            self.assertEqual(result['epochs_completed'], 3)
            self.assertEqual(result['best_epoch_zero_based'], 1)
            self.assertEqual(result['best_val_loss'], 2.5)

    def test_rejects_nonfinite_loss(self):
        with tempfile.TemporaryDirectory() as directory:
            history = Path(directory) / 'history.csv'
            output = Path(directory) / 'summary.json'
            history.write_text('epoch,loss,val_loss\n0,1.0,nan\n')
            with self.assertRaisesRegex(ValueError, 'Non-finite'):
                run(history, output)


if __name__ == '__main__':
    unittest.main()
