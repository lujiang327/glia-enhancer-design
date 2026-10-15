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
import pyBigWig
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from build_celltype_browser_tracks import main, count_endpoints, compare_mg_reference
from prepare_browser_tracks import digest


class CelltypeBrowser(unittest.TestCase):
    def test_four_donor_export_matches_existing_mg_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);inputs=root/'input';inputs.mkdir()
            sizes=root/'sizes';sizes.write_text('chr1\t130\n')
            donors={}
            for name in ['LGS1','LGS2','LGS3','LVG1']:
                path=inputs/(name+'.fragments.bedpe.gz')
                with gzip.open(path,'wt') as f:f.write('chr1\t10\t50\nchr1\t60\t120\n')
                donors[name]={'fragments':2,'cells':1,'sha256':digest(path)}
            reference=root/'raw.bw'
            with pyBigWig.open(str(reference),'w') as bw:
                bw.addHeader([('chr1',130)])
                bw.addEntries(['chr1']*4,[10,50,60,120],ends=[11,51,61,121],values=[4.]*4)
            config=root/'config.json'
            config.write_text(json.dumps({'bin_bp':100,'chrom_sizes_sha256':digest(sizes),
                'mg_reference_sha256':digest(reference),'groups':{'MG':{'keep_donor_tracks':True,'donors':donors}}}))
            output=root/'output'
            argv=['build','--config',str(config),'--group','MG','--pseudobulk-dir',str(inputs),
                '--chrom-sizes',str(sizes),'--output-dir',str(output),'--mg-reference',str(reference)]
            with patch.object(sys,'argv',argv),contextlib.redirect_stdout(io.StringIO()):main()
            report=json.loads((output/'browser_track_summary.json').read_text())
            self.assertEqual(report['fragments'],8)
            self.assertEqual(report['insertions'],16)
            self.assertEqual(report['mg_reference_check'],'PASS_EXACT_ALL_GENOME_INTEGER_BIN_MATCH')
            self.assertEqual(len(list(output.glob('*.bw'))),5)
            with pyBigWig.open(str(output/'MG_pooled.IPM_per_base.100bp.bw')) as bw:
                np.testing.assert_allclose(bw.values('chr1',0,1),[12/100*1e6/16],rtol=1e-6)
                np.testing.assert_allclose(bw.values('chr1',120,121),[4/30*1e6/16],rtol=1e-6)
            with self.assertRaises(FileExistsError),patch.object(sys,'argv',argv):main()

    def test_endpoint_convention_and_fragment_count_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'d.gz'
            with gzip.open(path,'wt') as f:f.write('chr1\t99\t100\nchr1\t1\t99\n')
            counts=count_endpoints(path,{'chr1':101},100,2)
            np.testing.assert_array_equal(counts['chr1'],[3,1])
            with self.assertRaises(ValueError):count_endpoints(path,{'chr1':100},100,2)
            with self.assertRaises(ValueError):count_endpoints(path,{'chr1':101},100,3)
            reference=Path(directory)/'different.bw'
            with pyBigWig.open(str(reference),'w') as bw:
                bw.addHeader([('chr1',101)])
                bw.addEntries(['chr1'],[1],ends=[2],values=[4.])
            with self.assertRaises(ValueError):compare_mg_reference(counts,reference,{'chr1':101},100)

if __name__=='__main__': unittest.main()
