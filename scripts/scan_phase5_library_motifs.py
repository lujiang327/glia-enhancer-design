#!/usr/bin/env python3
"""Descriptive both-strand FIMO scan; no TF-binding or motif causality claims."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
import pandas as pd


def annotate_hits(hits, library):
    lookup = library.set_index('sequence_id')
    if not set(hits.sequence_name).issubset(lookup.index):
        raise ValueError('Unknown FIMO sequence identifier')
    hits = hits.copy()
    hits['sequence_start_0based'] = hits.start.astype(int)-1
    hits['sequence_end_0based'] = hits.stop.astype(int)
    if ((hits.sequence_start_0based < 0) | (hits.sequence_end_0based > 500) |
        (hits.sequence_end_0based <= hits.sequence_start_0based)).any():
        raise ValueError('Invalid FIMO coordinates')
    metadata = lookup.loc[hits.sequence_name].reset_index(drop=True)
    for c in ['peak_id', 'sequence_role', 'chrom']:
        hits[c] = metadata[c].values
    hits['reference_start_0based'] = metadata.start.values + hits.sequence_start_0based.values
    hits['reference_end_0based'] = metadata.start.values + hits.sequence_end_0based.values
    positions = metadata.parent_position_0based.values
    hits['overlaps_edited_base'] = ((positions >= hits.sequence_start_0based.values) &
                                   (positions < hits.sequence_end_0based.values))
    return hits


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    a = p.parse_args()
    cfg = json.loads(a.config.read_text())
    for name in ['motifs', 'library', 'fasta']:
        if hashlib.sha256(Path(cfg[name]).read_bytes()).hexdigest() != cfg[name+'_sha256']:
            raise ValueError('Pinned input hash mismatch: '+name)
    a.output_dir.mkdir(parents=True, exist_ok=False)
    cmd = ['fimo', '--text', '--thresh', str(cfg['p_value_threshold']), cfg['motifs'], cfg['fasta']]
    with (a.output_dir/'fimo_hits.tsv').open('w') as out, (a.output_dir/'fimo.stderr.log').open('w') as err:
        subprocess.run(cmd, stdout=out, stderr=err, check=True)
    library = pd.read_csv(cfg['library'], sep='\t')
    hits = pd.read_csv(a.output_dir/'fimo_hits.tsv', sep='\t', comment='#')
    hits = annotate_hits(hits, library)
    hits.to_csv(a.output_dir/'annotated_motif_hits.tsv.gz', sep='\t', index=False)
    groups = {k: g for k, g in hits.groupby('sequence_name')}
    keycols = ['motif_id', 'sequence_start_0based', 'sequence_end_0based', 'strand']
    def sites(seq):
        h = groups.get(seq)
        return set() if h is None else set(h[keycols].itertuples(index=False, name=None))
    summaries = []
    for row in library.itertuples():
        h = groups.get(row.sequence_id)
        parent_sites = sites('parent__'+row.peak_id)
        current = sites(row.sequence_id)
        changed = current.symmetric_difference(parent_sites)
        if row.sequence_role != 'parent' and any(not (x[1] <= row.parent_position_0based < x[2]) for x in changed):
            raise ValueError('Motif change outside edited base')
        summaries.append(dict(sequence_id=row.sequence_id, motif_hit_count=len(current),
            motif_ids=';'.join(sorted(set() if h is None else set(h.motif_id))),
            motif_hits_overlapping_edit=0 if h is None else int(h.overlaps_edited_base.sum()),
            motif_sites_gained_at_threshold=len(current-parent_sites),
            motif_sites_lost_at_threshold=len(parent_sites-current)))
    result = library.drop(columns='motif_evidence_status').merge(pd.DataFrame(summaries), on='sequence_id', validate='one_to_one')
    result['motif_evidence_status'] = 'FIMO sequence match only; per-site attribution not assessed'
    result.to_csv(a.output_dir/'ranked_library_with_motif_matches.tsv.gz', sep='\t', index=False)
    summary = dict(status='MOTIF_SCAN_COMPLETE_SCIENTIFIC_REVIEW_PENDING', sequences=len(result), hits=len(hits),
                   command=cmd, fimo_version=subprocess.check_output(['fimo', '--version'], text=True).strip(),
                   config=cfg, threshold_scope='nominal site p-value; no q-values in text mode',
                   limitations=['matches do not establish TF binding or causal motif contribution',
                                'database includes enzyme motifs and redundant TF family matches; retain IDs for review',
                                'gained/lost means threshold crossing only; variant attribution not computed'])
    (a.output_dir/'motif_scan_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2), flush=True)

if __name__ == '__main__':
    main()
