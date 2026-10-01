from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from score_phase4_parent_sequences import reverse_complement, sha256_text


class Phase4ParentScoring(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ACGTN"), "NACGT")
        self.assertEqual(reverse_complement(reverse_complement("AACCGGTTN")), "AACCGGTTN")

    def test_sequence_hash_is_stable(self):
        self.assertEqual(
            sha256_text("ACGT"),
            "1dff3e84fe7877e0673b69bbddcf40124e396e3f9943dd890c91b6a09adb9af0",
        )


if __name__ == "__main__":
    unittest.main()
