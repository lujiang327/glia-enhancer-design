#!/usr/bin/env python3
"""Organize the complete draft into parent-centered groups, without selecting oligos."""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--library', type=Path, default=Path('reports/phase5/motif_attribution/ranked_library_with_parent_attribution.tsv.gz'))
    p.add_argument('--output-dir', type=Path, required=True)
    a = p.parse_args()
    library = pd.read_csv(a.library, sep='\t')
    parents = library[library.sequence_role.eq('parent')].copy()
    if len(parents) != 1000 or parents.peak_id.duplicated().any():
        raise ValueError('Expected 1,000 unique natural parents')
    cols = ['candidate_rank','peak_id','chrom','start','end','genomic_class','nearest_gene',
        'signed_distance_to_nearest_tss','parent_sequence','parent_sequence_sha256',
        'primary_logFC','primary_FDR','primary_donors_positive','primary_minimum_donor_logFC',
        'minimum_pairwise_logFC','minimum_glial_lodo_logFC','maximum_glial_lodo_FDR',
        'muller_mean_cpm','maximum_offtarget_mean_cpm','maximum_offtarget_cell_type','observed_specificity_log2_ratio',
        'prediction_forward_logcount','prediction_reverse_complement_logcount','prediction_mean_logcount',
        'parent_count_orientation_flag','attribution_qc_pass','attribution_qc_failure_reasons',
        'repeat_overlap_fraction','repeat_classes','repeat_families','umap_k50_parent_mean','umap_k100_parent_mean',
        'reference_umap_k50_mean_lt_0_8_review_flag','reference_umap_k100_mean_lt_0_8_review_flag',
        'repeat_overlap_ge_50pct_review_flag']
    grouped = parents[cols].rename(columns={c:'parent_'+c for c in cols if c.startswith('prediction_')})
    fields = ['sequence_id','sequence','sequence_sha256','parent_position_0based','genomic_position_0based','ref','alt',
        'prediction_forward_logcount','prediction_reverse_complement_logcount','prediction_mean_logcount',
        'delta_forward_logcount','delta_reverse_complement_logcount','delta_mean_logcount','robust_effect_lower_bound',
        'motif_sites_gained_at_threshold','motif_sites_lost_at_threshold',
        'counts_consensus_edited_parent_base_attribution','counts_consensus_edited_base_absolute_attribution_percentile',
        'counts_forward_alt_minus_ref_hypothetical_delta','counts_reverse_alt_minus_ref_hypothetical_delta',
        'edit_base_repeat_classes','edit_base_repeat_families','edit_base_reference_umap_k50','edit_base_reference_umap_k100']
    for role in ['gain','loss_control']:
        edits = library[library.sequence_role.eq(role)]
        if edits.peak_id.duplicated().any() or not set(edits.peak_id).issubset(set(parents.peak_id)):
            raise ValueError('Duplicate or orphan edits')
        grouped = grouped.merge(edits[['peak_id']+fields].rename(columns={c:role+'_'+c for c in fields}),
                                on='peak_id',how='left',validate='one_to_one')
        grouped['has_'+role] = grouped[role+'_sequence_id'].notna()
    grouped['unique_insert_count'] = 1 + grouped.has_gain.astype(int) + grouped.has_loss_control.astype(int)
    if grouped.unique_insert_count.sum() != len(library):
        raise ValueError('Library coverage mismatch')
    ids = set(parents.sequence_id)
    for role in ['gain','loss_control']:
        ids.update(grouped[role+'_sequence_id'].dropna())
    if ids != set(library.sequence_id):
        raise ValueError('Sequence identifier mismatch')
    for row in grouped.itertuples():
        for role in ['gain','loss_control']:
            if getattr(row,'has_'+role):
                seq = getattr(row,role+'_sequence')
                pos = int(getattr(row,role+'_parent_position_0based'))
                expected = row.parent_sequence[:pos]+getattr(row,role+'_alt')+row.parent_sequence[pos+1:]
                if seq != expected or hashlib.sha256(seq.encode()).hexdigest()!=getattr(row,role+'_sequence_sha256'):
                    raise ValueError('Parent/variant integrity mismatch')
    grouped['assembly'] = 'hg38'
    grouped['insert_length_bp'] = 500
    grouped['selection_status'] = 'complete draft; no final experimental subset selected'
    grouped['variant_offtarget_specificity_available'] = False
    grouped['synthesis_ready'] = False
    a.output_dir.mkdir(parents=True,exist_ok=False)
    grouped.sort_values('candidate_rank').to_csv(a.output_dir/'parent_gain_control_groups.tsv.gz',sep='\t',index=False)
    summary = dict(status='GROUPED_DRAFT_AWAITING_EXPERIMENTAL_PARAMETERS', parent_loci=len(grouped), unique_inserts=len(library),
        complete_triplets=int((grouped.has_gain & grouped.has_loss_control).sum()),
        parent_and_gain_only=int((grouped.has_gain & ~grouped.has_loss_control).sum()),
        parent_and_loss_control_only=int((~grouped.has_gain & grouped.has_loss_control).sum()),
        parent_only=int((~grouped.has_gain & ~grouped.has_loss_control).sum()),
        inserts_before_barcodes_and_adapters=len(library), required_inputs=['assay/vector','construct budget',
            'insert-length constraints','cloning/adaptor sequences and prohibited sites','barcode and replicate design'],
        source_sha256=hashlib.sha256(a.library.read_bytes()).hexdigest())
    (a.output_dir/'experimental_group_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    main()
