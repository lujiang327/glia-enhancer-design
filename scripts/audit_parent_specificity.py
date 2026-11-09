#!/usr/bin/env python3
"""Compare existing donor-aware DA statistics with the delivered browser tracks.

This reviews the fixed selected set. It does not refit tests or reselect parents.
"""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig

ROOT = Path(__file__).resolve().parents[1]
PAIRS = ['Astrocyte', 'Microglia', 'Rod', 'Cone', 'ON_cone_bipolar',
         'OFF_cone_bipolar', 'Rod_bipolar', 'GABA_amacrine', 'Gly_amacrine',
         'AII_amacrine', 'Horizontal', 'RGC']
DONORS = ['LGS1', 'LGS2', 'LGS3', 'LVG1']


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_track(path, parents):
    with pyBigWig.open(str(path)) as bw:
        # Exact values, not approximate zoom-level summaries. Keep partial bins.
        result = np.array([np.mean(bw.values(r.chrom, int(r.start), int(r.end)))
                           for r in parents.itertuples()])
    if not np.isfinite(result).all() or (result < 0).any():
        raise ValueError('Nonfinite or negative track signal: ' + str(path))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tracks-dir', type=Path, default=ROOT / 'data/intermediate/browser_tracks/retina_hg38_review_100bp')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'reports/browser_tracks/parent_specificity_audit')
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    source = ROOT / 'reports/phase3/candidate_ranking/top_1000_candidates.tsv.gz'
    parents = pd.read_csv(source, sep='\t').sort_values('candidate_rank').reset_index(drop=True)
    assert len(parents) == 1000 and parents.peak_id.is_unique
    assert parents.candidate_rank.tolist() == list(range(1, 1001))
    assert parents.end.sub(parents.start).eq(500).all()
    # Reconcile the candidate BED actually shipped to reviewers.
    bed = pd.read_csv(args.tracks_dir / 'MG_top1000_parents.hg38.bed', sep='\t', header=None)
    expected = {(r.chrom, r.start, r.end, f'rank{r.candidate_rank:04d}__{r.peak_id}') for r in parents.itertuples()}
    assert set(map(tuple, bed.iloc[:, :4].to_numpy())) == expected
    manifest = pd.read_csv(args.tracks_dir / 'track_index.tsv', sep='\t')
    pooled = pd.DataFrame()
    track_sources = []
    for group in ['MG'] + PAIRS:
        name = group + '_pooled.IPM_per_base.100bp.bw'
        path = args.tracks_dir / name
        actual = sha(path)
        assert actual == manifest.loc[manifest.file.eq(name), 'sha256'].item(), name
        track_sources.append(dict(file=name, sha256=actual))
        pooled[group] = read_track(path, parents)
    donors = pd.DataFrame()
    for donor in DONORS:
        values = {}
        for group in ['MG', 'Astrocyte', 'Microglia']:
            name = group + '_' + donor + '.IPM_per_base.100bp.bw'
            path = args.tracks_dir / name
            actual = sha(path)
            assert actual == manifest.loc[manifest.file.eq(name), 'sha256'].item(), name
            track_sources.append(dict(file=name, sha256=actual))
            values[group] = read_track(path, parents)
        for group, signal in values.items():
            donors[group + '_' + donor] = signal
        donors[donor + '_MG_above_both_glia'] = values['MG'] > np.maximum(values['Astrocyte'], values['Microglia'])
    off = pooled[PAIRS]
    best = off.max(axis=1)
    # No pseudocount needed for the strongest off-target among 12 here.
    assert best.gt(0).all()
    ratio = pooled.MG / best
    pair_logfc = parents[[ct + '_logFC' for ct in PAIRS]]
    pair_fdr = parents[[ct + '_FDR' for ct in PAIRS]]
    assert np.isfinite(pair_logfc).all().all()
    assert pair_fdr.ge(0).all().all() and pair_fdr.le(1).all().all()
    ids = parents[['candidate_rank', 'peak_id', 'chrom', 'start', 'end']].copy()
    ids['IGV_parent_locus'] = parents.chrom + ':' + (parents.start + 1).astype(str) + '-' + parents.end.astype(str)
    summary = parents[['candidate_rank', 'peak_id', 'genomic_class', 'nearest_gene',
                       'primary_logFC', 'primary_FDR', 'primary_donors_positive',
                       'primary_minimum_donor_logFC', 'maximum_glial_lodo_FDR']].copy()
    summary.insert(2, 'IGV_parent_locus', ids.IGV_parent_locus)
    summary['minimum_pairwise_log2FC'] = pair_logfc.min(axis=1)
    summary['maximum_pairwise_FDR'] = pair_fdr.max(axis=1)
    summary['significant_positive_pairwise_count_FDR_lt_0_05'] = ((pair_logfc.to_numpy() > 0) & (pair_fdr.to_numpy() < .05)).sum(axis=1)
    summary['browser_MG_mean_IPM_per_base'] = pooled.MG
    summary['browser_maximum_offtarget_mean_IPM_per_base'] = best
    summary['browser_strongest_offtarget'] = off.idxmax(axis=1)
    summary['browser_MG_to_strongest_offtarget_ratio'] = ratio
    summary['browser_donors_MG_above_both_glia_count'] = donors[[d + '_MG_above_both_glia' for d in DONORS]].sum(axis=1)
    summary.to_csv(out / 'parents_review.tsv', sep='\t', index=False)
    fc = pd.concat([ids, parents[['primary_logFC'] + [ct + '_logFC' for ct in PAIRS]]], axis=1)
    fdr = pd.concat([ids, parents[['primary_FDR'] + [ct + '_FDR' for ct in PAIRS]]], axis=1)
    signal = pd.concat([ids, pooled], axis=1)
    donor_export = pd.concat([ids, donors], axis=1)
    for name, frame in [('log2FC', fc), ('FDR', fdr), ('browser_signal', signal), ('donor_browser_signal', donor_export)]:
        frame.to_csv(out / (name + '.tsv'), sep='\t', index=False)
    contrasts = []
    for group, lc, fd in [('Rest (equal cell-type weights)', parents.primary_logFC, parents.primary_FDR)] + [(ct, parents[ct+'_logFC'], parents[ct+'_FDR']) for ct in PAIRS]:
        contrasts.append(dict(comparator=group, significant_positive_FDR_lt_0_05=int(((lc>0)&(fd<.05)).sum()),
                              minimum_log2FC=float(lc.min()), median_log2FC=float(lc.median()),
                              maximum_log2FC=float(lc.max()), median_FDR=float(fd.median()), maximum_FDR=float(fd.max())))
    contrast_frame = pd.DataFrame(contrasts)
    contrast_frame.to_csv(out / 'contrast_summary.tsv', sep='\t', index=False)
    weakest = int(summary.loc[ratio.idxmin(), 'candidate_rank'])
    review_ranks = [1, 500, 1000, weakest]
    result = dict(parents=1000, pairwise_comparisons=12,
        all_twelve_positive_and_FDR_lt_0_05=int(summary.significant_positive_pairwise_count_FDR_lt_0_05.eq(12).sum()),
        browser_MG_above_all_twelve_offtargets=int((pooled.MG > best).sum()),
        browser_ratio_min=float(ratio.min()), browser_ratio_median=float(ratio.median()), browser_ratio_max=float(ratio.max()),
        browser_ratio_ge_2=int(ratio.ge(2).sum()), browser_ratio_ge_4=int(ratio.ge(4).sum()),
        all_eight_glial_lodo_FDR_lt_0_05=int(parents.maximum_glial_lodo_FDR.lt(.05).sum()),
        all_four_browser_donors_MG_above_both_glia=int(summary.browser_donors_MG_above_both_glia_count.eq(4).sum()),
        per_donor_browser_MG_above_both_glia={d:int(donors[d+'_MG_above_both_glia'].sum()) for d in DONORS},
        weakest_browser_ratio_rank=weakest, browser_figure_ranks=review_ranks,
        track_sha256_verified=len(track_sources), candidate_BED_identity_check='PASS_ALL_1000',
        no_new_differential_tests=True, no_variant_specificity_evidence=True)
    (out/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    (out/'source_checksums.json').write_text(json.dumps(dict(candidate_table_sha256=sha(source), tracks=track_sources), indent=2)+'\n')
    # Compact scientific plots, separate from the spreadsheet's data tables.
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False})
    fig, axes = plt.subplots(1, 3, figsize=(16, 6.3))
    axes[0].boxplot([parents[ct+'_logFC'] for ct in PAIRS], vert=False, tick_labels=[ct.replace('_',' ') for ct in PAIRS], showfliers=False)
    axes[0].axvline(0, color='grey', ls='--'); axes[0].invert_yaxis()
    axes[0].set(xlabel='Donor-blocked MG / comparator log2FC', title='A  Effect sizes for all 1,000 parents')
    im=axes[1].imshow(-np.log10(pair_fdr.to_numpy()), aspect='auto', vmin=0, vmax=15, cmap='viridis')
    axes[1].set(xticks=range(12),xticklabels=[ct.replace('_',' ') for ct in PAIRS],ylabel='Original candidate rank',title='B  Existing pairwise significance')
    axes[1].tick_params(axis='x', rotation=90, labelsize=8)
    fig.colorbar(im, ax=axes[1], label='-log10(FDR), color capped at 15', shrink=.75)
    axes[2].hist(ratio, bins=30, color='#2878a5')
    axes[2].axvline(1,color='grey',ls='--'); axes[2].axvline(ratio.median(), color='black', ls=':')
    axes[2].set(xlabel='MG / strongest off-target mean browser signal',ylabel='Parent regions',title='C  Actual 100bp-binned BigWig values')
    fig.suptitle('Natural parents show MG-enriched accessibility with shared off-target signal',fontweight='bold',fontsize=16)
    fig.text(.03,.02,'Selected regions, not independent validation. FDRs use the original 147,359-peak universe per contrast. Browser ratios use exact 500bp parent means.\nBrowser tracks use pooled global insertion scaling; DA uses donor-blocked fragment counts and TMM. Accessibility enrichment does not establish reporter or variant specificity.',fontsize=9)
    fig.tight_layout(rect=[.01,.13,.99,.93])
    fig.savefig(out/'specificity_audit.png',dpi=180); plt.close(fig)
    # Include a deliberately weak case and fixed ranks, not only attractive peaks.
    fig, axes = plt.subplots(3, 4, figsize=(16, 7), sharex='col', sharey='col')
    groups=['MG','Astrocyte','Microglia']; colors=['#222222','#d88132','#2878a5']
    handles={g:pyBigWig.open(str(args.tracks_dir/(g+'_pooled.IPM_per_base.100bp.bw'))) for g in groups}
    try:
        for col, rank in enumerate(review_ranks):
            row=parents.iloc[rank-1]; start=int(row.start); end=int(row.end)
            left=start-1000; right=end+1000; center=(start+end)/2
            x=np.arange(left,right)-center
            window_max=0
            for i,(g,color) in enumerate(zip(groups,colors)):
                ax=axes[i,col]; values=np.asarray(handles[g].values(row.chrom,left,right))
                assert np.isfinite(values).all()
                window_max=max(window_max,float(values.max()))
                ax.fill_between(x,values,step='mid',color=color,alpha=.8)
                ax.axvspan(-250,250,color='#bbbbbb',alpha=.22)
                ax.set_ylabel(g if col==0 else '')
                ax.tick_params(labelsize=8)
                ax.set_xlim(-1250,1250)
            axes[0,col].set_title(f'Rank {rank}'+(' (weakest browser ratio)' if rank==weakest else '')+f'\n{ids.iloc[rank-1].IGV_parent_locus}\nMG / strongest off-target = {ratio.iloc[rank-1]:.2f}',fontsize=9)
            axes[0,col].set_ylim(0,window_max*1.12)
            axes[2,col].set_xlabel('Distance from parent center (bp)', fontsize=9)
        fig.suptitle('Same y-axis within each locus: a visible off-target peak can coexist with MG enrichment',fontweight='bold',fontsize=15)
        fig.text(.03,.025,'Shading is the exact 500bp parent. Each column has one common y-axis across MG, astrocyte and microglia; different loci may have different scales.\nObserved IPM per base, 100bp bins. Ratios compare MG with the strongest of all 12 off-targets. Includes prespecified ranks 1, 500, 1000 and the weakest ratio.',fontsize=9)
        fig.tight_layout(rect=[.01,.12,.99,.93],h_pad=.5,w_pad=2)
        fig.savefig(out/'shared_scale_locus_examples.png',dpi=180); plt.close(fig)
    finally:
        for handle in handles.values(): handle.close()
    # Intermediate JSON preserves numbers as numbers for the workbook author.
    tables={'parents':summary,'log2FC':fc,'FDR':fdr,'browser':signal,'contrasts':contrast_frame}
    # Avoid DataFrame.to_json, which rounds extremely small FDRs to zero.
    payload=dict(summary=result, comparisons=['Rest']+PAIRS, tables={})
    for name,frame in tables.items():
        payload['tables'][name]=dict(columns=list(frame.columns), rows=[
            [None if pd.isna(v) else v.item() if isinstance(v,np.generic) else v for v in row]
            for row in frame.itertuples(index=False,name=None)])
    (out/'workbook_inputs.json').write_text(json.dumps(payload,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
