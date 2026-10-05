#!/usr/bin/env python3
"""Verify motif attribution outputs and summarize descriptive parent context."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def contextual_flags(sites):
    result = sites.copy()
    result['enzyme_motif_reference_match'] = result.motif_id.str.match(r'^(DNASE|TN5)_')
    result['nonenzyme_positive_parent_count_context'] = (
        ~result.enzyme_motif_reference_match & result.counts_site_positive_both_orientations &
        result.counts_site_top10pct_parent_bases.gt(0))
    result['nonenzyme_negative_parent_count_context'] = (
        ~result.enzyme_motif_reference_match & result.counts_site_negative_both_orientations &
        result.counts_site_top10pct_parent_bases.gt(0))
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--results-dir', type=Path, required=True)
    a = p.parse_args()
    root = a.results_dir
    for line in (root/'output_checksums.sha256').read_text().splitlines():
        expected, name = line.split()
        assert hashlib.sha256((root/name).read_bytes()).hexdigest() == expected
    library = pd.read_csv(root/'ranked_library_with_parent_attribution.tsv.gz', sep='\t')
    sites = pd.read_csv(root/'edit_motif_parent_attribution.tsv.gz', sep='\t')
    original = pd.read_csv('reports/phase5/library_motifs/ranked_library_with_motif_matches.tsv.gz', sep='\t')
    old_sites = pd.read_csv('reports/phase5/library_motifs/edit_overlapping_motif_comparisons.tsv.gz', sep='\t')
    columns = [c for c in original if c != 'motif_evidence_status']
    pd.testing.assert_frame_equal(library[columns], original[columns], check_exact=False, rtol=1e-12, atol=1e-12)
    columns = [c for c in old_sites if c != 'attribution_overlap_status']
    pd.testing.assert_frame_equal(sites[columns], old_sites[columns], check_exact=False, rtol=1e-12, atol=1e-12)
    assert len(library) == library.sequence_id.nunique() == 2884
    assert len(sites) == 59828
    edits = library[library.sequence_role.ne('parent')].copy()
    parents = library[library.sequence_role.eq('parent')]
    for head in ['counts', 'profile']:
        for orientation in ['forward', 'reverse', 'consensus']:
            prefix = head+'_'+orientation+'_'
            assert parents[prefix+'alt_minus_ref_hypothetical_delta'].isna().all()
            assert np.isfinite(edits[prefix+'alt_minus_ref_hypothetical_delta']).all()
            assert edits[prefix+'edited_base_absolute_attribution_percentile'].between(0,1).all()
        assert np.isfinite(sites[head+'_site_sum_consensus']).all()
        assert sites[head+'_site_absolute_attribution_fraction_of_parent'].between(0,1).all()
        width = sites.sequence_end_0based - sites.sequence_start_0based
        assert (sites[head+'_site_top10pct_parent_bases'] <= np.minimum(width,50)).all()
        expected_positive = sites[head+'_site_sum_forward'].gt(0) & sites[head+'_site_sum_reverse_complement_aligned'].gt(0)
        expected_negative = sites[head+'_site_sum_forward'].lt(0) & sites[head+'_site_sum_reverse_complement_aligned'].lt(0)
        assert expected_positive.equals(sites[head+'_site_positive_both_orientations'])
        assert expected_negative.equals(sites[head+'_site_negative_both_orientations'])
    for orientation in ['forward', 'reverse']:
        direction = edits.sequence_role.map({'gain':1, 'loss_control':-1})
        assert (edits['counts_'+orientation+'_alt_minus_ref_hypothetical_delta']*direction > 0).all()
    sites = contextual_flags(sites)
    summary = dict(status='PASS_INTEGRITY_DESCRIPTIVE_PARENT_CONTEXT_REVIEW', source_job=63285123,
        elapsed='00:01:01', max_rss_kb=5326092, sequences=len(library), edits=len(edits), motif_comparisons=len(sites),
        enzyme_motif_comparisons=int(sites.enzyme_motif_reference_match.sum()), roles={},
        inference_scope='parent attribution context; no causal TF or variant specificity assignment',
        independent_validation=False,
        independent_validation_reason='edits were nominated using these same parent attribution maps')
    for role, g in edits.groupby('sequence_role'):
        ss = sites[sites.sequence_role.eq(role)]
        has_positive = set(ss[ss.nonenzyme_positive_parent_count_context].sequence_id)
        has_lost_positive = set(ss[ss.nonenzyme_positive_parent_count_context & ss.match_state_at_nominal_threshold.eq('lost')].sequence_id)
        has_lost_negative = set(ss[ss.nonenzyme_negative_parent_count_context & ss.match_state_at_nominal_threshold.eq('lost')].sequence_id)
        summary['roles'][role] = dict(sequences=len(g), count_delta_direction_agrees_both_orientations=int(g.counts_hypothetical_edit_direction_agrees_both_orientations.sum()),
            median_edited_parent_base_absolute_attribution_percentile=float(g.counts_consensus_edited_base_absolute_attribution_percentile.median()),
            median_direct_delta_mean_logcount=float(g.delta_mean_logcount.median()),
            median_robust_effect_lower_bound=float(g.robust_effect_lower_bound.median()),
            edits_with_nonenzyme_positive_parent_count_context=len(has_positive),
            edits_with_lost_match_in_nonenzyme_positive_parent_count_context=len(has_lost_positive),
            edits_with_lost_match_in_nonenzyme_negative_parent_count_context=len(has_lost_negative),
            hypothetical_vs_direct_effect_spearman=float(g.counts_consensus_alt_minus_ref_hypothetical_delta.corr(g.delta_mean_logcount,method='spearman')))
    sites.to_csv(root/'edit_motif_attribution_review.tsv.gz', sep='\t', index=False)
    # Review examples in original DA rank order, not a new functional ranking.
    example_ids = set(edits.sort_values('candidate_rank').groupby('sequence_role').head(10).sequence_id)
    examples = sites[sites.sequence_id.isin(example_ids)].sort_values(
        ['candidate_rank','sequence_role','enzyme_motif_reference_match','counts_site_absolute_attribution_fraction_of_parent'],
        ascending=[True,True,True,False]).groupby('sequence_id',sort=False).head(5)
    examples.to_csv(root/'example_edit_motif_context.tsv.gz', sep='\t', index=False)
    (root/'attribution_review_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    fig, axes = plt.subplots(1,2, figsize=(10,4), constrained_layout=True)
    for role, g in edits.groupby('sequence_role'):
        axes[0].hist(g.counts_consensus_edited_base_absolute_attribution_percentile, bins=np.linspace(0,1,26), alpha=.55, label=role)
        axes[1].scatter(g.counts_consensus_alt_minus_ref_hypothetical_delta,g.delta_mean_logcount,s=8,alpha=.35,label=role)
    axes[0].set(xlabel='Edited parent-base absolute count attribution percentile',ylabel='Number of selected edits')
    axes[1].set(xlabel='Parent hypothetical alternate-minus-reference count attribution',ylabel='Direct mean model logcount change')
    for ax in axes: ax.legend(frameon=False)
    fig.savefig(root/'motif_attribution_review.png', dpi=180)
    plt.close(fig)
    print(json.dumps(summary,indent=2))

if __name__ == '__main__':
    main()
