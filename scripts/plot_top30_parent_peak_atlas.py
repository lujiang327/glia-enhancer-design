#!/usr/bin/env python3
"""Visualize the most consistently significant 30 natural MG parent regions.

Selection uses existing DA statistics only, before reading/flanking track shapes.
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
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
PAIRS = ['Astrocyte','Microglia','Rod','Cone','ON_cone_bipolar',
         'OFF_cone_bipolar','Rod_bipolar','GABA_amacrine','Gly_amacrine',
         'AII_amacrine','Horizontal','RGC']
GROUPS = ['MG'] + PAIRS
LABELS = {'MG':'Müller glia', 'ON_cone_bipolar':'ON cone bipolar',
          'OFF_cone_bipolar':'OFF cone bipolar','Rod_bipolar':'Rod bipolar',
          'GABA_amacrine':'GABA amacrine','Gly_amacrine':'Gly amacrine',
          'AII_amacrine':'AII amacrine','RGC':'Ganglion cells'}
COLORS = {'MG':'#202020','Astrocyte':'#D88132','Microglia':'#2878A5'}


def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()


def activity_selection(parents, config, output_dir):
    """Activity-first ranking within explicit specificity and donor safeguards."""
    audit=ROOT/'reports/browser_tracks/parent_specificity_audit/parents_review.tsv'
    evidence=pd.read_csv(audit,sep='\t',float_precision='round_trip')
    fields=['candidate_rank','peak_id','browser_MG_mean_IPM_per_base',
            'browser_MG_to_strongest_offtarget_ratio','browser_donors_MG_above_both_glia_count']
    p=parents.merge(evidence[fields],on=['candidate_rank','peak_id'],validate='one_to_one')
    assert len(p)==1000
    e=config['eligibility']
    base=(p.worst_pairwise_FDR.lt(e['maximum_pairwise_FDR_lt']) &
          p.minimum_pairwise_logFC.ge(e['minimum_pairwise_log2FC_gte']) &
          p.maximum_glial_lodo_FDR.lt(e['maximum_glial_leave_one_donor_out_FDR_lt']) &
          p.primary_donors_positive.eq(e['primary_positive_donors']) &
          p.browser_donors_MG_above_both_glia_count.eq(e['browser_donors_MG_above_both_glia']))
    eligible=base & p.browser_MG_to_strongest_offtarget_ratio.ge(e['pooled_browser_MG_to_strongest_offtarget_gte'])
    order=['browser_MG_mean_IPM_per_base','browser_MG_to_strongest_offtarget_ratio','candidate_rank']
    ranked=p.sort_values(order,ascending=[False,False,True])
    ranked['activity_specificity_eligible']=eligible.loc[ranked.index]
    ranked.to_csv(output_dir/'all_1000_activity_review.tsv',sep='\t',index=False)
    selected=ranked[ranked.activity_specificity_eligible].head(config['target_count']).copy().reset_index(drop=True)
    assert len(selected)==30
    old=pd.read_csv(ROOT/'reports/browser_tracks/top30_peak_atlas/top30_selected_parents.tsv',sep='\t')
    old=p[p.peak_id.isin(old.peak_id)].copy()
    def metrics(g):
        return dict(regions=len(g),median_MG_browser_signal=float(g.browser_MG_mean_IPM_per_base.median()),
                    median_equal_donor_TMM_MG_CPM=float(g.muller_mean_cpm.median()),
                    minimum_browser_specificity=float(g.browser_MG_to_strongest_offtarget_ratio.min()),
                    median_browser_specificity=float(g.browser_MG_to_strongest_offtarget_ratio.median()))
    comparison={'old_significance_top30':metrics(old),'new_activity_specificity_top30':metrics(selected),
                'overlap_parents':len(set(old.peak_id)&set(selected.peak_id)),
                'eligible_parents':int(eligible.sum()),'original_universe':1000,
                'browser_audit_sha256':digest(audit)}
    comparison['median_MG_signal_change_fraction']=comparison['new_activity_specificity_top30']['median_MG_browser_signal']/comparison['old_significance_top30']['median_MG_browser_signal']-1
    activity_only=ranked.head(30)
    comparison['activity_only_top30_without_additional_filters']=metrics(activity_only)
    comparison['maximum_possible_top30_median_MG_signal_gain_in_fixed_1000']=float(activity_only.browser_MG_mean_IPM_per_base.median()/old.browser_MG_mean_IPM_per_base.median()-1)
    sensitivity=[]
    for threshold in config['threshold_sensitivity_browser_ratios']:
        keep=base & p.browser_MG_to_strongest_offtarget_ratio.ge(threshold)
        group=ranked.loc[keep.loc[ranked.index]].head(30)
        sensitivity.append(dict(minimum_browser_ratio=threshold,eligible=int(keep.sum()),**metrics(group)))
    pd.DataFrame(sensitivity).to_csv(output_dir/'threshold_sensitivity.csv',index=False)
    (output_dir/'ranking_comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
    fig,axes=plt.subplots(1,3,figsize=(15,5.8))
    ax=axes[0]
    ax.scatter(p.browser_MG_mean_IPM_per_base,p.browser_MG_to_strongest_offtarget_ratio,s=10,color='#BBC4CA',alpha=.5,label='All 1,000')
    ax.scatter(old.browser_MG_mean_IPM_per_base,old.browser_MG_to_strongest_offtarget_ratio,s=45,facecolors='none',edgecolors='#D88132',label='Original significance top 30')
    ax.scatter(selected.browser_MG_mean_IPM_per_base,selected.browser_MG_to_strongest_offtarget_ratio,s=22,color='#2878A5',label='New activity-first top 30')
    ax.axhline(3,color='#7C8791',ls='--',lw=1)
    ax.set(xlabel='Observed MG mean signal (IPM/base)',ylabel='MG / strongest off-target signal',title='A  Accessibility and specificity tradeoff')
    ax.legend(frameon=False,fontsize=8)
    for ax,metric,title,ylabel in [
        (axes[1],'browser_MG_mean_IPM_per_base','B  MG signal increases modestly','MG mean signal (IPM/base)'),
        (axes[2],'browser_MG_to_strongest_offtarget_ratio','C  Specificity is lower but remains >=3','MG / strongest off-target signal')]:
        ax.boxplot([old[metric],selected[metric]],tick_labels=['Significance\ntop 30','Activity-first\ntop 30'],showfliers=False)
        for j,g in enumerate([old,selected]):
            ax.scatter(np.full(len(g),j+1)+np.linspace(-.12,.12,len(g)),g[metric],s=13,color=['#D88132','#2878A5'][j],alpha=.7)
        ax.set(title=title,ylabel=ylabel)
    fig.suptitle('A stronger-activity review list, with an explicit specificity floor',weight='bold',fontsize=16)
    fig.text(.035,.025,'Same 1,000-parent universe and normalization. Selection metrics are not independent validation. Signal is pooled, depth-normalized ATAC, not reporter activity.\nNew review filters also require positive significant pairwise effects and donor consistency. A stricter 4-fold browser floor reduces the accessible activity range.',fontsize=9)
    fig.tight_layout(rect=[0,.13,1,.91]);fig.savefig(output_dir/'activity_vs_significance.png',dpi=190);plt.close(fig)
    return selected,comparison


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--tracks-dir',type=Path,default=ROOT/'data/intermediate/browser_tracks/retina_hg38_review_100bp')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'reports/browser_tracks/top30_peak_atlas')
    parser.add_argument('--pdf',type=Path,default=ROOT/'output/pdf/MG_top30_parent_peak_atlas.pdf')
    parser.add_argument('--ranking',choices=['significance','activity_specificity'],default='significance')
    parser.add_argument('--ranking-config',type=Path,default=ROOT/'config/activity_specificity_top30_review.json')
    a=parser.parse_args()
    active=a.ranking=='activity_specificity'
    if active and a.output_dir.resolve()==(ROOT/'reports/browser_tracks/top30_peak_atlas').resolve():
        parser.error('Activity review needs a separate --output-dir and --pdf to preserve the original atlas')
    if active and a.pdf.resolve()==(ROOT/'output/pdf/MG_top30_parent_peak_atlas.pdf').resolve():
        parser.error('Do not overwrite the original significance atlas PDF')
    a.output_dir.mkdir(parents=True,exist_ok=True)
    a.pdf.parent.mkdir(parents=True,exist_ok=True)
    source=ROOT/'reports/phase3/candidate_ranking/top_1000_candidates.tsv.gz'
    all_parents=pd.read_csv(source,sep='\t',float_precision='round_trip')
    assert len(all_parents)==1000 and all_parents.peak_id.is_unique
    fcols=[ct+'_FDR' for ct in PAIRS]; lcols=[ct+'_logFC' for ct in PAIRS]
    assert np.isfinite(all_parents[fcols+lcols]).all().all()
    assert all_parents[fcols].ge(0).all().all() and all_parents[fcols].le(1).all().all()
    all_parents['worst_pairwise_FDR']=all_parents[fcols].max(axis=1)
    all_parents['weakest_FDR_comparator']=all_parents[fcols].idxmax(axis=1).str.removesuffix('_FDR')
    comparison=None;ranking_config=None
    if active:
        ranking_config=json.loads(a.ranking_config.read_text())
        selected,comparison=activity_selection(all_parents,ranking_config,a.output_dir)
    else:
        selected=all_parents.sort_values(['worst_pairwise_FDR','candidate_rank']).head(30).copy().reset_index(drop=True)
    selected.insert(0,'display_order',np.arange(1,31))
    selected['IGV_parent_locus']=selected.chrom+':'+(selected.start+1).astype(str)+'-'+selected.end.astype(str)
    assert (selected[fcols]<.05).all().all() and (selected[lcols]>0).all().all()
    assert selected.end.sub(selected.start).eq(500).all()
    manifest=pd.read_csv(a.tracks_dir/'track_index.tsv',sep='\t')
    handles={}; source_hashes={}; means=np.zeros((30,13)); windows=[]
    try:
        for g in GROUPS:
            p=a.tracks_dir/(g+'_pooled.IPM_per_base.100bp.bw')
            source_hashes[p.name]=digest(p)
            assert source_hashes[p.name]==manifest.loc[manifest.file.eq(p.name),'sha256'].item()
            handles[g]=pyBigWig.open(str(p))
        for i,row in selected.iterrows():
            start,end=int(row.start),int(row.end)
            left,right=start-1000,end+1000
            assert left>=0 and right<=handles['MG'].chroms(row.chrom)
            values={g:np.asarray(b.values(row.chrom,left,right)) for g,b in handles.items()}
            assert all(np.isfinite(v).all() and (v>=0).all() for v in values.values())
            for j,g in enumerate(GROUPS):means[i,j]=values[g][1000:1500].mean()
            windows.append(dict(left=left,right=right,x=np.arange(left,right)-(start+end)/2,values=values))
    finally:
        for b in handles.values():b.close()
    assert (means[:,0]>means[:,1:].max(axis=1)).all()
    selected['browser_MG_to_strongest_offtarget']=means[:,0]/means[:,1:].max(axis=1)
    selected['browser_strongest_offtarget']=[PAIRS[i] for i in means[:,1:].argmax(axis=1)]
    if active:
        assert np.allclose(means[:,0],selected.browser_MG_mean_IPM_per_base,rtol=1e-13,atol=0)
        assert np.allclose(selected.browser_MG_to_strongest_offtarget,selected.browser_MG_to_strongest_offtarget_ratio,rtol=1e-13,atol=0)
        compact=selected[['display_order','candidate_rank','peak_id','IGV_parent_locus','nearest_gene','genomic_class',
            'browser_MG_mean_IPM_per_base','muller_mean_cpm','browser_strongest_offtarget',
            'browser_MG_to_strongest_offtarget','minimum_pairwise_logFC','worst_pairwise_FDR',
            'maximum_glial_lodo_FDR','browser_donors_MG_above_both_glia_count']].copy()
        compact.rename(columns={'display_order':'activity_rank','candidate_rank':'original_candidate_rank'}).to_csv(a.output_dir/'top30_activity_specificity.csv',index=False)
    selected.to_csv(a.output_dir/'top30_selected_parents.tsv',sep='\t',index=False)
    signal=selected[['display_order','candidate_rank','peak_id','IGV_parent_locus']].copy()
    for j,g in enumerate(GROUPS):signal[g+'_mean_IPM_per_base']=means[:,j]
    signal.to_csv(a.output_dir/'top30_browser_signal.tsv',sep='\t',index=False)
    # Save the actual displayed samples for exact plot reproduction without BigWigs.
    np.savez_compressed(a.output_dir/'displayed_track_values.npz',
                        values=np.stack([[w['values'][g] for g in GROUPS] for w in windows]),
                        relative_positions=windows[0]['x'], groups=np.asarray(GROUPS))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    figures=[]
    fig=plt.figure(figsize=(18,11.7))
    gs=fig.add_gridspec(1,3,left=.105,right=.95,bottom=.24,top=.845,width_ratios=[13,12,3],wspace=.23)
    axes=[fig.add_subplot(gs[0,i]) for i in range(3)]
    title='The 30 strongest MG ATAC parents passing specificity and donor safeguards' if active else 'The 30 parents with the strongest significance across all off-target comparisons'
    subtitle='Ranked by MG mean 500bp signal; require >=3-fold observed specificity, minimum fitted log2FC >=1.5, and donor safeguards.' if active else 'Selected by the smallest worst pairwise FDR among the 1,000 parents. Track shapes did not determine selection.'
    fig.text(.035,.955,title,fontsize=19,weight='bold')
    fig.text(.035,.918,subtitle,fontsize=11,color='#435361')
    rowlabels=[f'{i+1:02d}  parent #{int(r.candidate_rank)}' for i,r in selected.iterrows()]
    im=axes[0].imshow(means,aspect='auto',cmap='magma',vmin=0,vmax=means.max(),interpolation='nearest')
    axes[0].set(title='A  Observed accessibility in the 500bp parent',yticks=range(30),yticklabels=rowlabels,
                xticks=range(13),xticklabels=[LABELS.get(g,g) for g in GROUPS])
    axes[0].tick_params(axis='x',rotation=65,labelsize=9)
    for label in axes[0].get_xticklabels():label.set_ha('right')
    axes[0].tick_params(axis='y',labelsize=9)
    p=axes[0].get_position(); cax=fig.add_axes([p.x0,.095,p.width,.015])
    fig.colorbar(im,cax=cax,orientation='horizontal',label='Mean insertion density (IPM/base); one scale for all cells and regions')
    effects=selected[lcols].to_numpy()
    im2=axes[1].imshow(effects,aspect='auto',cmap='Blues',vmin=0,vmax=12,interpolation='nearest')
    axes[1].set(title='B  MG / off-target fitted log2 fold change',yticks=[],xticks=range(12),xticklabels=[LABELS.get(g,g) for g in PAIRS])
    axes[1].tick_params(axis='x',rotation=65,labelsize=9)
    for label in axes[1].get_xticklabels():label.set_ha('right')
    p=axes[1].get_position(); cax=fig.add_axes([p.x0,.095,p.width,.015])
    fig.colorbar(im2,cax=cax,orientation='horizontal',label='Donor-blocked log2FC; 1 = 2-fold, 2 = 4-fold, 3 = 8-fold')
    for i in range(30):
        for j in range(12):
            axes[1].text(j,i,f'{effects[i,j]:.1f}',ha='center',va='center',fontsize=7,color='white' if effects[i,j]>6.5 else '#24313C')
    axes[2].set_xlim(0,1); axes[2].set_ylim(29.5,-.5);axes[2].axis('off')
    axes[2].set_title('C  Worst FDR',fontsize=11)
    for i,r in selected.iterrows():axes[2].text(.5,i,f'{r.worst_pairwise_FDR:.1e}',ha='center',va='center',fontsize=9)
    fig.text(.035,.045,'Read across each row: brighter observed MG signal, positive pairwise effects, and small FDRs support enrichment. Shared glial signal is still present.',fontsize=10)
    fig.text(.035,.022,'Selected best cases, not representative of all 1,000 or independent validation. Natural-locus accessibility does not establish reporter or variant specificity.  |  1 / 7',fontsize=9,color='#56616B')
    overview=a.output_dir/'top30_visual_overview.png';fig.savefig(overview,dpi=200,facecolor='white');plt.close(fig);figures.append(overview)
    global_ymax=max(float(v.max()) for w in windows for v in w['values'].values())*1.10
    for page in range(6):
        fig=plt.figure(figsize=(18,11.7))
        gs=fig.add_gridspec(13,5,left=.10,right=.98,bottom=.15,top=.775,hspace=.20,wspace=.25)
        fig.text(.035,.956,f'Observed peak atlas: regions {page*5+1:02d}-{page*5+5:02d} of the selected 30',fontsize=20,weight='bold')
        scale_note='One common y-axis across ALL 30 regions and ALL 13 cell types. Gray shading marks the exact 500bp parent.' if active else 'One common y-axis across all 13 cell types within each column. Gray shading marks the exact 500bp parent.'
        legend_note='MG: black     Astrocytes: orange     Microglia: blue     Neurons: gray. Activity-first ranking; no independent autoscaling.' if active else 'MG: black     Astrocytes: orange     Microglia: blue     Neurons: gray. Columns may use different scales; compare down a column.'
        fig.text(.035,.918,scale_note,fontsize=11)
        fig.text(.035,.888,legend_note,fontsize=10,color='#435361')
        for c in range(5):
            i=page*5+c;row=selected.iloc[i];window=windows[i]
            ymax=global_ymax if active else max(float(v.max()) for v in window['values'].values())*1.10
            assert ymax>0
            for ri,g in enumerate(GROUPS):
                ax=fig.add_subplot(gs[ri,c])
                v=window['values'][g]
                ax.axvspan(-250,250,color='#DCE4E9',alpha=.6,zorder=0)
                ax.fill_between(window['x'],v,step='mid',color=COLORS.get(g,'#929BA2'),linewidth=.3,zorder=2)
                ax.set_xlim(-1250,1250);ax.set_ylim(0,ymax)
                ax.spines[['left','right','top']].set_visible(False)
                ax.spines['bottom'].set_color('#D6DDE2')
                ax.set_yticks([])
                if c==0:
                    ax.set_ylabel(LABELS.get(g,g),rotation=0,ha='right',va='center',labelpad=9,fontsize=9,weight='bold' if g=='MG' else 'normal',color=COLORS.get(g,'#435361'))
                if ri<12:ax.set_xticks([])
                else:
                    ax.set_xticks([-1000,0,1000]);ax.tick_params(axis='x',labelsize=8)
                    ax.set_xlabel('bp from parent center',fontsize=9)
                if ri==0:
                    gene=row.nearest_gene if pd.notna(row.nearest_gene) else 'unassigned'
                    ax.set_title(f'{i+1:02d}  Parent #{int(row.candidate_rank)}  ({gene})\n{row.IGV_parent_locus}\nWorst FDR {row.worst_pairwise_FDR:.1e}; min log2FC {row.minimum_pairwise_logFC:.2f}\nBrowser MG / max off-target {row.browser_MG_to_strongest_offtarget:.2f}\nAll track y-ranges: 0 to {ymax:.4f}',fontsize=8.5,pad=14,linespacing=1.4)
                    ax.spines['left'].set_visible(True)
                    ax.set_yticks([0,ymax]);ax.set_yticklabels(['0',f'{ymax:.3f}'],fontsize=7)
                    ax.tick_params(axis='y',length=2,pad=1)
        fig.text(.035,.082,'Signal: observed pooled fragment-end insertion density per million insertions per base, averaged over 100bp bins. No separate cell-type rescaling.',fontsize=10)
        fig.text(.035,.057,'Statistics: donor-blocked edgeR/TMM. Browser heights use global pooled scaling, so their ratios need not equal fitted fold changes. Four donors contributed.',fontsize=9,color='#56616B')
        fig.text(.035,.031,f'These are selected natural MG-enriched regions. Off-target peaks can remain. No reporter activity or designed-variant specificity is established.  |  {page+2} / 7',fontsize=9,color='#56616B')
        path=a.output_dir/f'top30_track_atlas_{page+1:02d}.png';fig.savefig(path,dpi=200,facecolor='white');plt.close(fig);figures.append(path)
    # Assemble the verified chart images into a portable, consistently sized PDF.
    pdf=canvas.Canvas(str(a.pdf),pagesize=(18*inch,11.7*inch))
    pdf.setTitle('Top 30 MG parent regions: visual accessibility and differential evidence')
    pdf.setAuthor('Glia regulatory-element project')
    for figure in figures:
        pdf.drawImage(str(figure),0,0,width=18*inch,height=11.7*inch)
        pdf.showPage()
    pdf.save()
    reader=PdfReader(str(a.pdf));assert len(reader.pages)==7
    for page in reader.pages:
        assert len(page.images)==1
        assert float(page.mediabox.width)==1296
    summary=dict(selection='Activity descending after specificity and donor filters' if active else '30 smallest maximum FDR across all 12 original pairwise contrasts; tie break original candidate rank',
        selection_used_track_shapes=False,parents=30,cell_types=13,pdf_pages=7,
        original_candidate_ranks=selected.candidate_rank.astype(int).tolist(),
        worst_pairwise_FDR_min=float(selected.worst_pairwise_FDR.min()),worst_pairwise_FDR_max=float(selected.worst_pairwise_FDR.max()),
        minimum_pairwise_log2FC=float(selected.minimum_pairwise_logFC.min()),
        browser_MG_to_strongest_offtarget_min=float(selected.browser_MG_to_strongest_offtarget.min()),
        browser_MG_to_strongest_offtarget_median=float(selected.browser_MG_to_strongest_offtarget.median()),
        browser_MG_to_strongest_offtarget_max=float(selected.browser_MG_to_strongest_offtarget.max()),
        all_glial_lodo_FDR_lt_0_05=int(selected.maximum_glial_lodo_FDR.lt(.05).sum()),
        normalization='Browser: global IPM/base 100bp means. DA: donor-blocked fragment overlaps with TMM.',
        source_table_sha256=digest(source),source_track_sha256=source_hashes,
        view='2.5kb per locus, centered on exact 500bp parent; all 13 y-ranges identical within each locus',
        scope='Selected best cases. Not representative of the entire library, independent validation, reporter evidence or variant specificity.')
    if active:
        summary.update(ranking_config=ranking_config,ranking_config_sha256=digest(a.ranking_config),ranking_comparison=comparison,
                       view='2.5kb per locus; one global y-range across every locus and cell type',atlas_global_ymax=global_ymax)
    (a.output_dir/'selection_and_provenance.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='source_track_sha256'},indent=2))


if __name__=='__main__':main()
