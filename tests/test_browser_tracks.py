from pathlib import Path
import sys
import unittest
import tempfile
import numpy as np
import pyBigWig
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from prepare_browser_tracks import bin_interval_counts, export


class BrowserBins(unittest.TestCase):
    def test_real_bigwig_export_conserves_normalized_signal(self):
        with tempfile.TemporaryDirectory() as directory:
            source=Path(directory)/'source.bw'
            destination=Path(directory)/'normalized.bw'
            with pyBigWig.open(str(source),'w') as bw:
                bw.addHeader([('chr1',250)])
                bw.addEntries(['chr1']*3,[10,99,249],ends=[11,100,250],values=[3.,2.,1.])
            result=export(source,destination,6,100)
            self.assertEqual(result['source_insertions'],6)
            with pyBigWig.open(str(destination)) as bw:
                np.testing.assert_allclose(bw.values('chr1',10,11),[5/100*1e6/6],rtol=1e-6)
                self.assertEqual(bw.values('chr1',150,151),[0.])
                np.testing.assert_allclose(bw.values('chr1',249,250),[1/50*1e6/6],rtol=1e-6)

    def test_spanning_intervals_and_partial_bin(self):
        bins=bin_interval_counts([(90,220,2),(200,250,1)],0,250,100)
        np.testing.assert_array_equal(bins,[20,200,90])
        self.assertEqual(bins.sum(),310)

    def test_chunk_boundary_is_not_double_counted(self):
        intervals=[(90,220,2),(200,250,1)]
        first=bin_interval_counts(intervals,0,100,100)
        second=bin_interval_counts(intervals,100,250,100)
        np.testing.assert_array_equal(np.concatenate([first,second]),[20,200,90])

    def test_empty_and_single_base_insertions(self):
        np.testing.assert_array_equal(bin_interval_counts(None,0,101,100),[0,0])
        np.testing.assert_array_equal(bin_interval_counts([(99,100,3),(100,101,2)],0,101,100),[3,2])

if __name__=='__main__': unittest.main()
