import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pyBigWig
from vendor.shift_math import get_ref_pwms,compute_shift_ATAC
from build_training_tracks import run as build_tracks
from prepare_training_regions import choose_bucket,sequence_ok


class PreparationContracts(unittest.TestCase):
    def test_shift_recovers_known_offset(self):
        root=Path(__file__).resolve().parents[1];rp,rm=get_ref_pwms(root/'config/ATAC.ref.motifs.txt')
        for startplus,startminus,expected in [(10,10,(4,-5)),(14,5,(0,0))]:
            p=np.ones((40,4))/4;m=p.copy()
            p[startplus:startplus+20]=np.mean(list(rp.values()),axis=0)
            m[startminus:startminus+20]=np.mean(list(rm.values()),axis=0)
            self.assertEqual(compute_shift_ATAC(rp,rm,p,m),expected)

    def test_duplicate_coordinates_across_donors_preserve_molecules(self):
        old=Path.cwd()
        with tempfile.TemporaryDirectory() as tmp:
            try:
                os.chdir(tmp)
                for d in ['reports/preparation','reports/qc','data/raw/reference','data/intermediate/donor_pseudobulks']:Path(d).mkdir(parents=True,exist_ok=True)
                Path('reports/preparation/fragment_shift.json').write_text(json.dumps(dict(status='PASS',plus_endpoint_delta=0,minus_endpoint_delta=0)))
                Path('data/raw/reference/hg38.canonical.chrom.sizes').write_text('chr1\t100\nchr2\t50\n')
                Path('reports/qc/donor_summary.tsv').write_text('donor\tretained_fragment_records\nA\t3\nB\t1\n')
                Path('data/intermediate/donor_pseudobulks/A.bedpe').write_text('chr1\t10\t20\nchr2\t1\t2\nchr1\t10\t21\n')
                Path('data/intermediate/donor_pseudobulks/B.bedpe').write_text('chr1\t10\t20\n')
                build_tracks()
                with pyBigWig.open('data/intermediate/training/muller.insertions.bw') as bw:
                    self.assertEqual(bw.values('chr1',10,11),[3.0])
                    self.assertEqual(bw.values('chr1',20,22),[2.0,1.0])
                    self.assertEqual(bw.header()['sumData'],8)
                    self.assertEqual(bw.values('chr2',1,3),[1.0,1.0])
            finally:os.chdir(old)

    def test_gc_fallback_is_bounded_and_context_checks(self):
        self.assertEqual(choose_bucket({45,50,60},49),50)
        with self.assertRaises(ValueError):choose_bucket(set(),49)
        self.assertFalse(sequence_ok('ACGTNACGT',0,5))
        self.assertFalse(sequence_ok('ACGT',-1,3))
        self.assertTrue(sequence_ok('ACGT',0,4))


if __name__=='__main__':unittest.main()
