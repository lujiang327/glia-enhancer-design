import json
from pathlib import Path
import tempfile
import unittest

from scripts.select_bias_motif_regions import run


class BiasMotifRegionTest(unittest.TestCase):
    def test_selection_is_reproducible_and_test_only(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            regions = root / 'regions.bed'
            fold = root / 'fold.json'
            output1 = root / 'one.bed'
            output2 = root / 'two.bed'
            manifest1 = root / 'one.json'
            manifest2 = root / 'two.json'
            regions.write_text(''.join(
                f'{chrom}\t{i}\t{i + 10}\tp{i}\t0\t.\t0\t0\t0\t5\n'
                for chrom in ('chr1', 'chr2') for i in range(10)
            ))
            fold.write_text(json.dumps({'train': ['chr2'], 'valid': [], 'test': ['chr1']}))
            run(regions, fold, 'test', 5, 42, output1, manifest1)
            run(regions, fold, 'test', 5, 42, output2, manifest2)
            self.assertEqual(output1.read_text(), output2.read_text())
            self.assertTrue(all(line.startswith('chr1\t') for line in output1.read_text().splitlines()))
            self.assertEqual(json.loads(manifest1.read_text())['selected_regions'], 5)


if __name__ == '__main__':
    unittest.main()
