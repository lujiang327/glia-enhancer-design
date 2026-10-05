#!/usr/bin/env python3
"""Assemble a review library; measured DA always describes the natural locus."""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd


def assemble(parents, edits, attribution):
    if parents.peak_id.duplicated().any() or edits.variant_id.duplicated().any():
        raise ValueError('Duplicate parent or variant identifiers')
    if attribution.peak_id.duplicated().any():
        raise ValueError('Duplicate attribution identifiers')
    parents = parents.merge(attribution, on=['peak_id', 'candidate_rank'], how='left',
                            validate='one_to_one', suffixes=('', '_attribution'), indicator=True)
    if not parents['_merge'].eq('both').all():
        raise ValueError('Missing attribution QC')
    parents = parents.drop(columns='_merge')
    lookup = parents.set_index('peak_id')
    records = []
    for _, row in parents.iterrows():
        d = row.to_dict()
        d.update(sequence_id='parent__' + row.peak_id, sequence_role='parent',
                 sequence=row.parent_sequence, sequence_sha256=row.parent_sequence_sha256,
                 prediction_forward_logcount=row.chrombpnet_forward_logcount,
                 prediction_reverse_complement_logcount=row.chrombpnet_reverse_complement_logcount,
                 prediction_mean_logcount=row.chrombpnet_mean_logcount,
                 delta_forward_logcount=0., delta_reverse_complement_logcount=0., delta_mean_logcount=0.)
        records.append(d)
    for _, edit in edits.iterrows():
        if edit.peak_id not in lookup.index:
            raise ValueError('Orphan variant')
        parent = lookup.loc[edit.peak_id]
        for c in ['candidate_rank', 'chrom', 'start', 'end', 'parent_sequence', 'parent_sequence_sha256']:
            if parent[c] != edit[c]:
                raise ValueError('Parent identity mismatch: ' + c)
        pos = int(edit.parent_position_0based)
        seq = parent.parent_sequence
        if not 0 <= pos < len(seq) or seq[pos] != edit.ref or edit.alt == edit.ref:
            raise ValueError('Invalid single-base substitution')
        if edit.variant_sequence != seq[:pos] + edit.alt + seq[pos+1:]:
            raise ValueError('Variant sequence mismatch')
        if int(edit.genomic_position_0based) != int(parent.start) + pos:
            raise ValueError('Edit coordinate mismatch')
        if not bool(edit.edit_screen_pass):
            raise ValueError('Unscreened edit')
        d = parent.to_dict()
        d.update(edit.to_dict())
        d.update(sequence_id=edit.variant_id, sequence_role=edit.design_class,
                 sequence=edit.variant_sequence, sequence_sha256=edit.variant_sequence_sha256,
                 prediction_forward_logcount=edit.variant_forward_logcount,
                 prediction_reverse_complement_logcount=edit.variant_reverse_complement_logcount,
                 prediction_mean_logcount=edit.variant_mean_logcount)
        records.append(d)
    result = pd.DataFrame(records)
    if result.sequence_id.duplicated().any():
        raise ValueError('Duplicate library identifiers')
    for row in result.itertuples():
        if len(row.sequence) != 500 or hashlib.sha256(row.sequence.encode()).hexdigest() != row.sequence_sha256:
            raise ValueError('Sequence length/hash mismatch')
    result['assembly'] = 'hg38'
    result['coordinate_convention'] = '0-based half-open; edited base 0-based'
    result['differential_accessibility_scope'] = 'natural parent locus only'
    result['prediction_scope'] = '500bp insert in its native 2114bp context; frozen MG model'
    result['variant_off_target_specificity_prediction_available'] = False
    result['predicted_variant_specificity_change'] = float('nan')
    result['motif_evidence_status'] = 'candidate sequence scan pending'
    result['library_status'] = 'DRAFT_REVIEW_PENDING_NOT_SYNTHESIS_READY'
    order = {'parent': 0, 'gain': 1, 'loss_control': 2}
    result['_order'] = result.sequence_role.map(order)
    return result.sort_values(['candidate_rank', '_order']).drop(columns='_order').reset_index(drop=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output-dir', type=Path, required=True)
    args = p.parse_args()
    root = Path('reports/phase4')
    paths = [root/'locus_annotation/parent_annotation_review.tsv.gz',
             root/'locus_annotation/edit_annotation_review.tsv.gz',
             root/'attribution/attribution_qc.tsv.gz']
    result = assemble(*(pd.read_csv(x, sep='\t') for x in paths))
    args.output_dir.mkdir(parents=True, exist_ok=False)
    result.to_csv(args.output_dir/'ranked_library_draft.tsv.gz', sep='\t', index=False)
    with (args.output_dir/'library_sequences.fa').open('w') as f:
        for row in result.itertuples():
            f.write('>' + row.sequence_id + '\n' + row.sequence + '\n')
    summary = dict(status='DRAFT_REVIEW_PENDING_NOT_SYNTHESIS_READY', sequences=len(result),
                   roles=result.sequence_role.value_counts().to_dict(), parents=result.peak_id.nunique(),
                   input_sha256={str(x): hashlib.sha256(x.read_bytes()).hexdigest() for x in paths},
                   pending=['candidate motif scan and review', 'off-target sequence predictions unavailable',
                            'final experimental selection and oligo/barcode/adapter design'])
    (args.output_dir/'library_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
