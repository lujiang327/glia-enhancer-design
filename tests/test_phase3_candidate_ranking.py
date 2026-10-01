import gzip
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from rank_phase3_candidates import annotate_peaks, percentiles, read_primary


class Phase3CandidateRanking(unittest.TestCase):
    def test_primary_parser_accepts_repeated_edger_coordinate_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "primary.tsv.gz"
            header = [
                "peak_id", "chrom", "start", "end", "peak_id", "chrom", "start", "end",
                "logFC", "logCPM", "F", "PValue", "FDR", "donor_logFC_LGS1",
                "donor_logFC_LGS2", "donor_logFC_LGS3", "donor_logFC_LVG1",
                "donors_positive", "minimum_donor_logFC",
            ]
            row = [
                "peak1", "chr1", "100", "600", "peak1", "chr1", "100", "600",
                "2.5", "3", "10", "0.001", "0.01", "1", "2", "3", "4", "4", "1",
            ]
            with gzip.open(path, "wt") as handle:
                handle.write("\t".join(header) + "\n" + "\t".join(row) + "\n")
            result = read_primary(path)
            self.assertEqual(result[0]["peak_id"], "peak1")
            self.assertEqual(result[0]["start"], 100)
            self.assertEqual(result[0]["primary_donors_positive"], 4)

    def test_percentiles_average_ties(self):
        self.assertEqual(percentiles([1, 2, 2, 4]), [0.0, 0.5, 0.5, 1.0])

    def test_annotation_priority(self):
        with tempfile.TemporaryDirectory() as tmp:
            refgene = Path(tmp) / "refgene.tsv.gz"
            rows = [
                "0\tNM_1\tchr1\t+\t1000\t5000\t1100\t4900\t2\t1000,4000,\t1200,5000,\t0\tGENE1\tcmpl\tcmpl\t0,0,",
            ]
            with gzip.open(refgene, "wt") as handle:
                handle.write("\n".join(rows) + "\n")
            peaks = [
                {"chrom": "chr1", "start": 900, "end": 1100},
                {"chrom": "chr1", "start": 4100, "end": 4200},
                {"chrom": "chr1", "start": 3000, "end": 3100},
                {"chrom": "chr1", "start": 8000, "end": 8100},
            ]
            result = annotate_peaks(peaks, refgene, promoter_window=500)
            self.assertEqual([item[0] for item in result], ["promoter", "exonic", "intronic", "distal_intergenic"])


if __name__ == "__main__":
    unittest.main()
