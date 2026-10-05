import unittest
from pathlib import Path
import sys
import tempfile
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from assemble_phase5_library import assemble
from scan_phase5_library_motifs import annotate_hits, read_fimo_hits


class LibraryContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[1]/'reports/phase4'
        cls.parents = pd.read_csv(root/'locus_annotation/parent_annotation_review.tsv.gz', sep='\t')
        cls.edits = pd.read_csv(root/'locus_annotation/edit_annotation_review.tsv.gz', sep='\t')
        cls.qc = pd.read_csv(root/'attribution/attribution_qc.tsv.gz', sep='\t')

    def test_all_natural_parents_and_da_evidence_retained(self):
        result = assemble(self.parents, self.edits, self.qc)
        self.assertEqual(len(result), 2884)
        natural = result[result.sequence_role.eq('parent')].set_index('peak_id')
        for row in self.parents.itertuples():
            self.assertEqual(natural.loc[row.peak_id, 'sequence'], row.parent_sequence)
        for row in result.itertuples():
            self.assertEqual(row.primary_FDR, natural.loc[row.peak_id, 'primary_FDR'])
        self.assertFalse(result.variant_off_target_specificity_prediction_available.any())
        self.assertTrue(result.predicted_variant_specificity_change.isna().all())

    def test_invalid_coordinate_or_sequence_rejected(self):
        for column in ['genomic_position_0based', 'variant_sequence']:
            edits = self.edits.copy()
            edits.loc[0, column] = edits.loc[0, column]+1 if column.startswith('genomic') else 'A'*500
            with self.assertRaises(ValueError):
                assemble(self.parents, edits, self.qc)

    def test_fimo_coordinates_and_edit_overlap_half_open(self):
        library = pd.DataFrame([dict(sequence_id='v', peak_id='p', sequence_role='gain', chrom='chr1', start=100,
                                     parent_position_0based=9)])
        hits = pd.DataFrame(dict(sequence_name=['v', 'v'], start=[1,11], stop=[10,20]))
        result = annotate_hits(hits, library)
        self.assertEqual(result.sequence_start_0based.tolist(), [0,10])
        self.assertEqual(result.reference_end_0based.tolist(), [110,120])
        self.assertEqual(result.overlaps_edited_base.tolist(), [True,False])
        hits.loc[0, 'stop'] = 501
        with self.assertRaises(ValueError):
            annotate_hits(hits, library)

    def test_legacy_and_modern_fimo_headers_preserve_first_match(self):
        headers = ['#pattern name\tsequence name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched sequence',
                   'motif_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence']
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'fimo.tsv'
            for header in headers:
                path.write_text(header+'\nNFI\tv\t1\t3\t+\t8.1\t1e-5\t\tACT\n')
                hits = read_fimo_hits(path)
                self.assertEqual(len(hits), 1)
                self.assertEqual(hits.iloc[0].motif_id, 'NFI')
                self.assertEqual(hits.iloc[0].sequence_name, 'v')
                path.write_text(header+'\n')
                self.assertTrue(read_fimo_hits(path).empty)
            path.write_text('unexpected\tcolumns\n')
            with self.assertRaises(ValueError):
                read_fimo_hits(path)

if __name__ == '__main__':
    unittest.main()
