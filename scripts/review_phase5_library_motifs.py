#!/usr/bin/env python3
"""Review recovered candidate motif matches without assigning causal TF identity."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--results-dir', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    a = p.parse_args()
    original = pd.read_csv('reports/phase5/library_draft/ranked_library_draft.tsv.gz', sep='\t')
    library = pd.read_csv(a.results_dir/'ranked_library_with_motif_matches.tsv.gz', sep='\t')
    original = original.set_index('sequence_id').loc[library.sequence_id].reset_index()
    columns = [c for c in original.columns if c != 'motif_evidence_status']
    pd.testing.assert_frame_equal(library[columns], original[columns], check_exact=False, rtol=1e-12, atol=1e-12)
    hits = pd.read_csv(a.results_dir/'annotated_motif_hits.tsv.gz', sep='\t')
    # FIMO rounds printed p-values; internally sub-threshold sites can print 1e-4.
    assert (hits['p-value'] <= 1e-4).all()
    assert hits['q-value'].isna().all()
    keys = ['motif_id', 'sequence_start_0based', 'sequence_end_0based', 'strand']
    assert not hits.duplicated(['sequence_name'] + keys).any()
    groups = {name: g.set_index(keys) for name, g in hits.groupby('sequence_name')}
    comparisons = []
    for row in library[library.sequence_role.ne('parent')].itertuples():
        parent = groups['parent__'+row.peak_id]
        variant = groups[row.sequence_id]
        pos = int(row.parent_position_0based)
        parent = parent[(parent.index.get_level_values(1) <= pos) & (parent.index.get_level_values(2) > pos)]
        variant = variant[(variant.index.get_level_values(1) <= pos) & (variant.index.get_level_values(2) > pos)]
        gained = lost = 0
        for key in sorted(set(parent.index) | set(variant.index)):
            ph = parent.loc[key] if key in parent.index else None
            vh = variant.loc[key] if key in variant.index else None
            state = 'retained' if ph is not None and vh is not None else ('gained' if vh is not None else 'lost')
            gained += state == 'gained'
            lost += state == 'lost'
            comparisons.append(dict(sequence_id=row.sequence_id, peak_id=row.peak_id,
                candidate_rank=row.candidate_rank, sequence_role=row.sequence_role,
                motif_id=key[0], sequence_start_0based=key[1], sequence_end_0based=key[2], strand=key[3],
                chrom=row.chrom, reference_start_0based=row.start+key[1], reference_end_0based=row.start+key[2],
                edit_position_0based=pos, match_state_at_nominal_threshold=state,
                parent_site_score=np.nan if ph is None else ph.score,
                variant_site_score=np.nan if vh is None else vh.score,
                parent_site_pvalue=np.nan if ph is None else ph['p-value'],
                variant_site_pvalue=np.nan if vh is None else vh['p-value'],
                delta_model_forward_logcount=row.delta_forward_logcount,
                delta_model_reverse_complement_logcount=row.delta_reverse_complement_logcount,
                attribution_overlap_status='pending; match is not causal attribution evidence'))
        assert gained == row.motif_sites_gained_at_threshold
        assert lost == row.motif_sites_lost_at_threshold
    comparisons = pd.DataFrame(comparisons)
    summary = dict(status='PASS_LIBRARY_AND_MOTIF_MATCH_INTEGRITY_ATTRIBUTION_OVERLAP_PENDING',
        sequences=len(library), hits=len(hits), edit_overlapping_motif_comparisons=len(comparisons), roles={})
    for name, g in library.groupby('sequence_role'):
        summary['roles'][name] = dict(sequences=len(g), median_matches=float(g.motif_hit_count.median()),
            sequences_with_gained_sites=int(g.motif_sites_gained_at_threshold.gt(0).sum()),
            sequences_with_lost_sites=int(g.motif_sites_lost_at_threshold.gt(0).sum()),
            sequences_with_any_threshold_change=int((g.motif_sites_gained_at_threshold+g.motif_sites_lost_at_threshold).gt(0).sum()))
    a.output_dir.mkdir(parents=True, exist_ok=False)
    library.to_csv(a.output_dir/'ranked_library_with_motif_matches.tsv.gz', sep='\t', index=False)
    comparisons.to_csv(a.output_dir/'edit_overlapping_motif_comparisons.tsv.gz', sep='\t', index=False)
    (a.output_dir/'motif_review_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    (a.output_dir/'motif_scan_summary.json').write_bytes((a.results_dir/'motif_scan_summary.json').read_bytes())
    for name in ['motif_scan_summary.json', 'annotated_motif_hits.tsv.gz', 'ranked_library_with_motif_matches.tsv.gz']:
        print(name, hashlib.sha256((a.results_dir/name).read_bytes()).hexdigest())
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
