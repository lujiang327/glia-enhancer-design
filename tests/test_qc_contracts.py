import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_fragments import audit
from pool_qc import merged_bed, overlaps
from model_gate import require_validated


class QCContracts(unittest.TestCase):
    def test_barcode_collisions_and_read_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);frag=root/'f.gz';rna=root/'r.gz'
            with gzip.open(frag,'wt') as f:f.write('chr1\t10\t20\tAAA-1\t100\nchr1\t30\t40\tAAA-1\t1\nchrM\t10\t20\tBBB-1\t2\n')
            with gzip.open(rna,'wt') as f:f.write('1_AAA-1\n2_AAA-1\n1_BBB-1\n')
            result=audit(frag,rna,root/'out')
            self.assertEqual(result['fragment_records'],3)
            self.assertEqual(result['read_support_sum'],103)
            self.assertEqual(result['records_by_number_of_matching_RNA_samples']['2'],2)
            self.assertEqual(result['chr1_22_X_records'],2)
            self.assertIsNone(result['biological_cell_count'])

    def test_blacklist_half_open_and_nested_intervals(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'b.bed';p.write_text('chr1\t10\t30\nchr1\t15\t20\nchr1\t40\t50\n')
            bed=merged_bed(p)
            self.assertFalse(overlaps(bed,'chr1',0,10))
            self.assertTrue(overlaps(bed,'chr1',25,31))
            self.assertFalse(overlaps(bed,'chr1',30,40))
            self.assertTrue(overlaps(bed,'chr1',49,51))

    def test_model_gate_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'g.json';p.write_text('{"status":"NOT_TRAINED"}')
            with self.assertRaises(ValueError):require_validated(p)
            p.write_text('{"status":"VALIDATED","checks":{}}')
            with self.assertRaises(ValueError):require_validated(p)

    def test_chromosome_splits(self):
        root=Path(__file__).resolve().parents[1];universe={f'chr{i}' for i in range(1,23)}|{'chrX'};tests=[]
        for p in sorted((root/'config/splits').glob('*.json')):
            d=json.loads(p.read_text());a,b,c=[set(d[k]) for k in ['train','valid','test']]
            self.assertFalse(a&b or a&c or b&c)
            self.assertEqual(a|b|c,universe);tests.extend(c)
        self.assertEqual(len(tests),len(universe));self.assertEqual(set(tests),universe)


if __name__=='__main__':unittest.main()
