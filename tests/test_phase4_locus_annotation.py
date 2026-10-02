import contextlib
import gzip
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import pyBigWig

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from annotate_phase4_loci import main, mappability_features, repeat_features, union_bp
from download_phase4_annotation_assets import digest


class LocusAnnotation(unittest.TestCase):
    def test_repeat_union_clips_and_does_not_count_boundary_touch(self):
        hits = [(50, 120, "r1", "SINE", "Alu"), (90, 160, "r2", "SINE", "Alu"),
                (200, 230, "r3", "LINE", "L1")]
        features = repeat_features(hits, 100, 200)
        self.assertEqual(features["repeat_overlap_bp"], 60)
        self.assertEqual(features["repeat_hit_count"], 2)
        self.assertEqual(features["repeat_classes"], "SINE")
        self.assertEqual(union_bp([(1, 5), (2, 4), (5, 8)]), 7)

    def test_omitted_umap_bases_contribute_zero(self):
        summary, values = mappability_features([1, np.nan, .5, np.nan])
        self.assertEqual(summary["mean"], .375)
        self.assertEqual(summary["zero_fraction"], .5)
        self.assertEqual(summary["omitted_fraction"], .5)
        self.assertEqual(values.tolist(), [1, 0, .5, 0])
        with self.assertRaises(ValueError):
            mappability_features([np.inf])

    def test_end_to_end_snapshot_repeat_and_bigwig_annotation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            parents = pd.DataFrame({"peak_id": ["p{}".format(i) for i in range(1000)],
                                    "chrom": "chr1", "start": np.arange(1000) * 600,
                                    "end": np.arange(1000) * 600 + 500,
                                    "parent_sequence": "A" * 500})
            parents.to_csv(root / "parents.tsv", sep="\t", index=False)
            edit = parents.iloc[:1].copy()
            edit["variant_id"] = "p0:gain:0A>C"
            edit["parent_position_0based"] = 0
            edit["genomic_position_0based"] = 0
            edit.to_csv(root / "edits.tsv", sep="\t", index=False)
            (root / "reference.fai").write_text("chr1\t600000\t0\t0\t0\n")
            columns = ["bin", "swScore", "milliDiv", "milliDel", "milliIns", "genoName",
                       "genoStart", "genoEnd", "genoLeft", "strand", "repName", "repClass",
                       "repFamily", "repStart", "repEnd", "repLeft", "id"]
            (root / "rmsk.sql").write_text("CREATE TABLE `rmsk` (\n" +
                "\n".join("  `{}` int,".format(name) for name in columns) + "\n);")
            with gzip.open(root / "rmsk.txt.gz", "wt") as handle:
                for start, end in ((0, 100), (50, 150), (500, 550), (590, 650)):
                    values = ["0", "1", "0", "0", "0", "chr1", str(start), str(end),
                              "0", "+", "AluTest", "SINE", "Alu", "0", "0", "0", "1"]
                    handle.write("\t".join(values) + "\n")
            for k, end, score in ((50, 50, .5), (100, 500, 1.0)):
                with pyBigWig.open(str(root / "k{}.Umap.MultiTrackMappability.bw".format(k)), "w") as bw:
                    bw.addHeader([("chr1", 600000)])
                    bw.addEntries(["chr1"], [0], ends=[end], values=[score])
            names = ["rmsk.sql", "rmsk.txt.gz", "k50.Umap.MultiTrackMappability.bw", "k100.Umap.MultiTrackMappability.bw"]
            manifest = {"assets": {name: {"url": "test:" + name, "sha256": digest(root / name)} for name in names}}
            (root / "asset_snapshot.json").write_text(json.dumps(manifest))
            config = {"assets_directory": str(root), "parents": str(root / "parents.tsv"),
                      "edits": str(root / "edits.tsv"), "reference_fai": str(root / "reference.fai"),
                      "assets": [{"name": name, "url": "test:" + name} for name in names],
                      "mappability": {"read_lengths_bp": [50, 100]}, "next_gate": "REVIEW"}
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config))
            with patch.object(sys, "argv", ["annotate", "--config", str(config_path), "--output-dir", str(root / "output")]), contextlib.redirect_stdout(io.StringIO()):
                main()
            result = pd.read_csv(root / "output/annotated_natural_parents.tsv.gz", sep="\t")
            self.assertEqual(len(result), 1000)
            self.assertEqual(result.iloc[0].repeat_overlap_bp, 150)
            self.assertEqual(result.iloc[1].repeat_overlap_bp, 50)
            self.assertAlmostEqual(result.iloc[0].umap_k50_parent_mean, .05)
            self.assertEqual(result.iloc[0].umap_k100_parent_mean, 1)
            self.assertEqual(result.iloc[1].umap_k100_parent_mean, 0)
            annotated_edits = pd.read_csv(root / "output/annotated_best_single_edits.tsv.gz", sep="\t")
            self.assertEqual(annotated_edits.iloc[0].edit_base_repeat_overlap_bp, 1)
            self.assertFalse(annotated_edits.iloc[0].variant_mappability_prediction_available)
            (root / "rmsk.sql").write_text("changed source")
            with patch.object(sys, "argv", ["annotate", "--config", str(config_path), "--output-dir", str(root / "second_output")]):
                with self.assertRaisesRegex(ValueError, "checksum or URL mismatch"):
                    main()


if __name__ == "__main__":
    unittest.main()
