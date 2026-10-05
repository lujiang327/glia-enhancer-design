#!/usr/bin/env python3
"""Connect motif windows to parent attribution; no variant SHAP or causal TF claims."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

BASES = 'ACGT'


def project(scores, one_hot):
    scores = np.asarray(scores, dtype=np.float64)
    one_hot = np.asarray(one_hot)
    if scores.shape != one_hot.shape or scores.ndim != 2 or scores.shape[1] != 4:
        raise ValueError('Expected matching length-by-ACGT arrays')
    if not np.isfinite(scores).all() or not np.isin(one_hot, [0, 1]).all() or not (one_hot.sum(axis=1) == 1).all():
        raise ValueError('Non-finite scores or ambiguous parent sequence')
    return (scores * one_hot).sum(axis=1)


def site_features(forward, reverse, consensus, start, end):
    if not 0 <= start < end <= len(consensus):
        raise ValueError('Motif interval outside parent')
    f = float(forward[start:end].sum())
    r = float(reverse[start:end].sum())
    total = float(np.abs(consensus).sum())
    top = np.argsort(-np.abs(consensus), kind='stable')[:max(1, int(np.ceil(len(consensus)*.1)))]
    return dict(site_sum_forward=f, site_sum_reverse_complement_aligned=r,
        site_sum_consensus=float(consensus[start:end].sum()),
        site_absolute_attribution_fraction_of_parent=float(np.abs(consensus[start:end]).sum()/total) if total else np.nan,
        site_top10pct_parent_bases=int(((top >= start) & (top < end)).sum()),
        site_positive_both_orientations=f > 0 and r > 0,
        site_negative_both_orientations=f < 0 and r < 0)


def edit_features(hypothetical, projected, position, ref, alt):
    if not 0 <= position < len(projected) or ref not in BASES or alt not in BASES or ref == alt:
        raise ValueError('Invalid edit')
    return dict(edited_parent_base_attribution=float(projected[position]),
        edited_base_absolute_attribution_percentile=float((np.abs(projected) <= abs(projected[position])).mean()),
        alt_minus_ref_hypothetical_delta=float(hypothetical[position, BASES.index(alt)]-hypothetical[position, BASES.index(ref)]))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--attributions', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    a = p.parse_args()
    cfg = json.loads(a.config.read_text())
    for name in ['library', 'motif_comparisons']:
        if hashlib.sha256(Path(cfg[name]).read_bytes()).hexdigest() != cfg[name+'_sha256']:
            raise ValueError('Input checksum mismatch: '+name)
    if hashlib.sha256(a.attributions.read_bytes()).hexdigest() != cfg['attribution_sha256']:
        raise ValueError('Attribution checksum mismatch')
    import h5py
    library = pd.read_csv(cfg['library'], sep='\t')
    comparisons = pd.read_csv(cfg['motif_comparisons'], sep='\t')
    if library.sequence_id.duplicated().any():
        raise ValueError('Duplicate sequence IDs')
    parents = library[library.sequence_role.eq('parent')].set_index('peak_id')
    edits = library[library.sequence_role.ne('parent')].set_index('sequence_id')
    cache = {}
    with h5py.File(a.attributions, 'r') as h:
        begin, finish = cfg['parent_slice_zero_based_half_open']
        if (int(h.attrs['parent_start']), int(h.attrs['parent_end'])) != (begin, finish):
            raise ValueError('Parent slice mismatch')
        ids = h['peak_id'].asstr()[:]
        ranks = h['candidate_rank'][:]
        if len(set(ids)) != len(ids) or set(ids) != set(parents.index):
            raise ValueError('Parent IDs mismatch')
        for i, pid in enumerate(ids):
            row = parents.loc[pid]
            one_hot = h['sequence_one_hot'][i, begin:finish]
            seq = ''.join(BASES[x] for x in one_hot.argmax(axis=1))
            if seq != row.parent_sequence or int(ranks[i]) != int(row.candidate_rank):
                raise ValueError('Parent sequence/rank mismatch')
            cache[pid] = {}
            for head in ['counts', 'profile']:
                maps = {}
                for orientation, dataset in [('forward', 'forward_hypothetical'),
                    ('reverse', 'reverse_complement_aligned_hypothetical'), ('consensus', 'consensus_hypothetical')]:
                    hypothetical = h[head+'/'+dataset][i, begin:finish].astype(np.float64)
                    maps[orientation] = (hypothetical, project(hypothetical, one_hot))
                cache[pid][head] = maps
    edit_rows = []
    for sid, row in edits.iterrows():
        pos = int(row.parent_position_0based)
        if row.parent_sequence[pos] != row.ref:
            raise ValueError('Edited reference base mismatch')
        features = {'sequence_id': sid}
        for head in ['counts', 'profile']:
            for orientation in ['forward', 'reverse', 'consensus']:
                hyp, projected = cache[row.peak_id][head][orientation]
                features.update({head+'_'+orientation+'_'+k: v for k, v in edit_features(hyp, projected, pos, row.ref, row.alt).items()})
            direction = 1 if row.sequence_role == 'gain' else -1
            features[head+'_hypothetical_edit_direction_agrees_both_orientations'] = all(
                direction * features[head+'_'+orientation+'_alt_minus_ref_hypothetical_delta'] > 0
                for orientation in ['forward', 'reverse'])
        edit_rows.append(features)
    site_rows = []
    for row in comparisons.itertuples():
        edit = edits.loc[row.sequence_id]
        if row.peak_id != edit.peak_id or not row.sequence_start_0based <= edit.parent_position_0based < row.sequence_end_0based:
            raise ValueError('Motif/edit identity or overlap mismatch')
        features = {}
        for head in ['counts', 'profile']:
            maps = cache[row.peak_id][head]
            features.update({head+'_'+k: v for k, v in site_features(maps['forward'][1], maps['reverse'][1],
                maps['consensus'][1], int(row.sequence_start_0based), int(row.sequence_end_0based)).items()})
        site_rows.append(features)
    sites = pd.concat([comparisons.reset_index(drop=True), pd.DataFrame(site_rows)], axis=1)
    sites['attribution_overlap_status'] = 'parent count/profile projected attribution in matched interval; no variant attribution'
    result = library.merge(pd.DataFrame(edit_rows), on='sequence_id', how='left', validate='one_to_one')
    result['motif_evidence_status'] = 'sequence matches with parent attribution context; no causal TF assignment'
    a.output_dir.mkdir(parents=True, exist_ok=False)
    sites.to_csv(a.output_dir/'edit_motif_parent_attribution.tsv.gz', sep='\t', index=False)
    result.to_csv(a.output_dir/'ranked_library_with_parent_attribution.tsv.gz', sep='\t', index=False)
    summary = dict(status='PARENT_ATTRIBUTION_OVERLAP_COMPLETE_SCIENTIFIC_REVIEW_PENDING', parents=len(parents),
        edits=len(edits), motif_comparisons=len(sites), config=cfg,
        counts_direction_agrees_both_orientations=int(pd.DataFrame(edit_rows).counts_hypothetical_edit_direction_agrees_both_orientations.sum()),
        limitations=['all attribution describes natural parent, including gained motif intervals',
            'float16 saved scores introduce quantization; exact edit effects remain direct float32 model predictions',
            'overlapping/redundant motifs share bases; their attribution sums must not be added together',
            'site sign and magnitude are descriptive, not statistical significance or causal TF evidence',
            'profile attribution concerns profile shape, not total accessibility direction'])
    (a.output_dir/'motif_attribution_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2), flush=True)

if __name__ == '__main__':
    main()
