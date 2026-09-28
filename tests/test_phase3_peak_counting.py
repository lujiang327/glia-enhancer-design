import gzip
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from count_phase3_consensus_peaks import count_fragments, load_peaks


class Phase3PeakCounting(unittest.TestCase):
    def test_half_open_overlap_and_multipeak_tracking(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            peak_path = root / "peaks.bed"
            peak_path.write_text(
                "chr1\t10\t20\tretina_peak_000001\n"
                "chr1\t20\t30\tretina_peak_000002\n"
                "chr2\t5\t15\tretina_peak_000003\n"
            )
            fragment_path = root / "fragments.gz"
            with gzip.open(str(fragment_path), "wt") as handle:
                handle.write("chr1\t0\t10\n")
                handle.write("chr1\t10\t11\n")
                handle.write("chr1\t19\t21\n")
                handle.write("chr1\t30\t40\n")
                handle.write("chr2\t14\t16\n")
            peaks, ordered = load_peaks(peak_path)
            counts, flow = count_fragments(fragment_path, peaks, len(ordered))
            self.assertEqual(counts, [2, 1, 1])
            self.assertEqual(flow["fragment_records"], 5)
            self.assertEqual(flow["fragments_in_any_peak"], 3)
            self.assertEqual(flow["peak_overlap_assignments"], 4)
            self.assertEqual(flow["fragments_overlapping_multiple_peaks"], 1)
            self.assertEqual(flow["extra_peak_assignments"], 1)

    def test_unsorted_fragments_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            peak_path = root / "peaks.bed"
            peak_path.write_text("chr1\t10\t20\tretina_peak_000001\n")
            fragment_path = root / "fragments.gz"
            with gzip.open(str(fragment_path), "wt") as handle:
                handle.write("chr1\t15\t16\nchr1\t14\t17\n")
            peaks, ordered = load_peaks(peak_path)
            with self.assertRaises(ValueError):
                count_fragments(fragment_path, peaks, len(ordered))


if __name__ == "__main__":
    unittest.main()
