#!/usr/bin/env python3
"""Build an auditable pre-experiment figure review from existing measured outputs."""
import base64
import hashlib
import io
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'reports/decision_review'
SOURCES = []
BLUE, ORANGE, GREY, RED = '#2878a5', '#d88132', '#7c8791', '#ba4545'


def load_json(name):
    SOURCES.append(name)
    return json.loads((ROOT/name).read_text())


def table(name):
    SOURCES.append(name)
    text = (ROOT/name).read_text() if not name.endswith('.gz') else None
    # A few legacy exports contain literal escaped delimiters.
    if text and '\\t' in text.splitlines()[0]:
        return pd.read_csv(io.StringIO(text.replace('\\t','\t').replace('\\r\\n','\n')), sep='\t')
    return pd.read_csv(ROOT/name, sep='\t')


def finish(fig, name, title, footer):
    fig.suptitle(title, fontsize=16, weight='bold', y=.99)
    fig.text(.05,.015,footer,fontsize=9,color='#414b54',va='bottom')
    fig.tight_layout(rect=[.025,.095,.99,.94],h_pad=3,w_pad=3)
    for extension in ['png','svg']:
        fig.savefig(OUT/(name+'.'+extension),dpi=180,facecolor='white')
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    seeds = pd.DataFrame(load_json('reports/training/chrombpnet_fold0_neg1_seed_reproducibility.json')['per_seed'])
    donor = load_json('reports/training/chrombpnet_fold0_seed42_neg1_donor_qc.json')['donors']
    bias = load_json('reports/training/bias_fold0_seed42_numerical_qc.json')['metrics']
    library = table('reports/phase5/motif_attribution/ranked_library_with_parent_attribution.tsv.gz')
    parents = library[library.sequence_role.eq('parent')].sort_values('candidate_rank')
    edits = library[library.sequence_role.ne('parent')]
    scan = load_json('reports/phase4/perturbation_full/perturbation_full_summary.json')
    review = load_json('reports/phase5/motif_attribution/attribution_review_summary.json')
    pseudo = table('reports/phase3/pseudobulk_qc.tsv')
    assert len(parents)==1000 and len(edits)==1884
    seed_export = seeds[['seed','peak_counts_pearsonr','peak_median_jsd','nonpeak_counts_pearsonr']].copy()
    seed_export['bias_peak_counts_pearsonr'] = bias['peak_counts_pearsonr']
    seed_export['bias_peak_median_jsd'] = bias['peak_median_jsd']
    seed_export['bias_nonpeak_counts_pearsonr'] = bias['nonpeak_counts_pearsonr']
    seed_export.to_csv(OUT/'model_seed_evidence.csv',index=False)
    donor_export = pd.DataFrame([dict(donor=k,
        peak_model_pearson=v['model_metrics']['peaks']['counts_pearsonr'],
        peak_bias_pearson=v['matched_scaled_bias_metrics']['peaks']['counts_pearsonr']) for k,v in donor.items()])
    donor_export.to_csv(OUT/'model_donor_evidence.csv',index=False)
    fig,axs=plt.subplots(2,2,figsize=(12,8))
    x=np.arange(len(seeds))
    for ax,metric,baseline,ylabel,title in [
        (axs[0,0],'peak_counts_pearsonr',bias['peak_counts_pearsonr'],'Count Pearson r (higher is better)','A  Peak prediction reproduces across three seeds'),
        (axs[0,1],'peak_median_jsd',bias['peak_median_jsd'],'Median profile JSD (lower is better)','B  Peak profiles improve over enzyme bias'),
        (axs[1,1],'nonpeak_counts_pearsonr',bias['nonpeak_counts_pearsonr'],'Count Pearson r (higher is better)','D  Background calibration remains a weakness')]:
        ax.scatter(x,seeds[metric],s=75,color=BLUE,zorder=3,label='Full model')
        ax.axhline(baseline,color=RED,ls='--',label='Bias-only baseline')
        for i,v in enumerate(seeds[metric]): ax.annotate(f'{v:.3f}',(i,v),xytext=(0,8),textcoords='offset points',ha='center')
        ax.set(xticks=x,xticklabels=seeds.seed.astype(str),xlabel='Training seed',ylabel=ylabel,title=title,ylim=(0,1))
        ax.legend(frameon=False,loc='upper right')
    ax=axs[1,0]; x=np.arange(len(donor_export))
    ax.bar(x-.18,donor_export.peak_bias_pearson,width=.36,color=GREY,label='Bias only')
    ax.bar(x+.18,donor_export.peak_model_pearson,width=.36,color=BLUE,label='Full model')
    ax.set(xticks=x,xticklabels=donor_export.donor,ylim=(0,1),ylabel='Peak count Pearson r',title='C  All four contributing donors improve over bias')
    ax.legend(frameon=False)
    finish(fig,'01_model_validation','1. The MG model learns useful sequence signal, with a background caveat',
        'All numerical tests shown use held-out chr1, chr3, chr6. Donors all contributed to training: this is donor consistency, not unseen-donor validation.\nSeed stability is aggregate performance; the designed edit effects have not yet been reproduced across seeds.')

    pairs=['Astrocyte','Microglia','Rod','Cone','ON_cone_bipolar','OFF_cone_bipolar','Rod_bipolar','GABA_amacrine','Gly_amacrine','AII_amacrine','Horizontal','RGC']
    fig,axs=plt.subplots(2,2,figsize=(12,9))
    ax=axs[0,0]
    for name,color in [('Microglia',BLUE),('Astrocyte',ORANGE),('Other',GREY)]:
        mask=parents.maximum_offtarget_cell_type.eq(name) if name!='Other' else ~parents.maximum_offtarget_cell_type.isin(['Microglia','Astrocyte'])
        ax.scatter(parents.loc[mask,'maximum_offtarget_mean_cpm'],parents.loc[mask,'muller_mean_cpm'],s=10,alpha=.55,color=color,label=f'{name} (n={mask.sum()})')
    ax.plot([.5,100],[.5,100],ls='--',color='black',lw=1,label='Equal accessibility')
    ax.set(xscale='log',yscale='log',xlim=(.5,100),ylim=(.5,100),xlabel='Highest off-target mean ATAC CPM',ylabel='MG mean ATAC CPM',title='A  Parents have higher measured MG accessibility')
    ax.legend(frameon=False,fontsize=8)
    ax=axs[0,1]
    image=ax.imshow(parents[[p+'_logFC' for p in pairs]].to_numpy(),aspect='auto',cmap='viridis',vmin=0,vmax=8,interpolation='nearest')
    ax.set(xticks=np.arange(len(pairs)),xticklabels=[p.replace('_',' ') for p in pairs],ylabel='Parent regions in original DA rank',title='B  Positive MG effect against all 12 off-targets')
    ax.tick_params(axis='x',rotation=65,labelsize=8)
    fig.colorbar(image,ax=ax,label='Donor-blocked MG vs cell-type log2FC')
    ax=axs[1,0]
    ax.hist(parents.primary_minimum_donor_logFC,bins=25,color=BLUE)
    ax.axvline(0,color='black',ls='--')
    ax.set(xlabel='Smallest MG-vs-rest log2FC among four donors',ylabel='Parent regions',title='C  Direction is consistent in all four donors')
    ax=axs[1,1]
    worst=np.sort(parents.maximum_glial_lodo_FDR.to_numpy())
    ax.plot(worst,np.arange(1,len(worst)+1)/len(worst),color=BLUE,lw=2)
    ax.axvline(.05,color=RED,ls='--')
    n=int((worst<.05).sum())
    ax.text(.05,.8,f'{n}/1,000 pass FDR < 0.05\nin every glial donor-omission refit',transform=ax.transAxes,fontsize=10)
    ax.set(xscale='log',xlabel='Worst FDR across eight glial donor-omission fits',ylabel='Fraction of selected parents',title='D  Donor sensitivity is favorable, with exceptions',ylim=(0,1.03))
    finish(fig,'02_parent_specificity','2. Natural parent regions have donor-supported MG-enriched accessibility',
        'These are selected regions: positive directions and specificity helped define eligibility/ranking. This is not independent experimental validation.\nATAC accessibility and statistical enrichment do not establish enhancer activity or mutant specificity.')
    parents[['peak_id','candidate_rank','muller_mean_cpm','maximum_offtarget_mean_cpm','maximum_offtarget_cell_type',
        'primary_minimum_donor_logFC','maximum_glial_lodo_FDR']+[p+'_logFC' for p in pairs]].to_csv(OUT/'parent_specificity_evidence.csv',index=False)

    fig,axs=plt.subplots(2,2,figsize=(12,8.5))
    ax=axs[0,0]
    categories=['Nominated','Both orientations\nintended effect','Effect/profile\nscreen passed','Best edit\nper parent']
    funnel=[]
    for role,color,offset in [('gain',BLUE,-.18),('loss_control',ORANGE,.18)]:
        data=scan['by_design_class'][role]
        values=[data['variants'],data['intended_effect_both_orientations'],data['screen_pass_variants'],int(edits.sequence_role.eq(role).sum())]
        ax.bar(np.arange(4)+offset,values,.36,color=color,label=role.replace('_',' '))
        funnel.extend(dict(role=role,stage=stage,edits=n) for stage,n in zip(categories,values))
        for i,n in enumerate(values): ax.text(i+offset,n+55,str(n),ha='center',fontsize=8)
    ax.set(xticks=range(4),xticklabels=categories,ylabel='Edits',title='A  Screening removes many nominated gain edits',ylim=(0,5500))
    ax.tick_params(axis='x',labelsize=8)
    ax.legend(frameon=False)
    for role,color in [('gain',BLUE),('loss_control',ORANGE)]:
        g=edits[edits.sequence_role.eq(role)]
        axs[0,1].scatter(g.delta_forward_logcount,g.delta_reverse_complement_logcount,s=8,alpha=.4,color=color,label=f'{role.replace("_"," ")} (n={len(g)})')
        axs[1,0].hist(g.counts_consensus_edited_base_absolute_attribution_percentile,bins=np.linspace(0,1,26),alpha=.55,color=color,label=role.replace('_',' '))
        axs[1,1].scatter(g.counts_consensus_alt_minus_ref_hypothetical_delta,g.delta_mean_logcount,s=8,alpha=.35,color=color)
    ax=axs[0,1];ax.axhline(0,color=GREY,lw=.8);ax.axvline(0,color=GREY,lw=.8)
    ax.set(xlabel='Forward direct logcount change',ylabel='Reverse-complement direct logcount change',title='B  Retained effects agree in both orientations')
    ax.legend(frameon=False,fontsize=8)
    axs[1,0].set(xlabel='Edited-base absolute parent attribution percentile',ylabel='Retained edits',title='C  Loss controls hit highly attributed parent bases')
    axs[1,0].legend(frameon=False)
    axs[1,1].set(xlabel='Hypothetical count-attribution change (alt - ref)',ylabel='Direct mean model logcount change',title='D  Gain-effect magnitude is less well explained')
    axs[1,1].text(.03,.04,'Within-class Spearman: gains 0.297; loss controls 0.665',transform=axs[1,1].transAxes,fontsize=9)
    finish(fig,'03_variant_predictions','3. Designs have screened model effects; their biological effects remain predictions',
        '947 attribution-eligible parents; 9,470 nominated substitutions. Retained variants were selected using these same attribution/model results.\nOrientation agreement is a selection criterion, not independent validation. Predictions use a native 2,114bp context, not an isolated 500bp reporter.')
    pd.DataFrame(funnel).to_csv(OUT/'variant_screen_evidence.csv',index=False)
    edits[['sequence_id','sequence_role','peak_id','delta_forward_logcount','delta_reverse_complement_logcount',
        'counts_consensus_edited_base_absolute_attribution_percentile','counts_consensus_alt_minus_ref_hypothetical_delta','delta_mean_logcount']].to_csv(OUT/'variant_effect_evidence.csv',index=False)

    fig,axs=plt.subplots(2,2,figsize=(12,8.5))
    ax=axs[0,0]
    largest=parents.maximum_offtarget_cell_type.value_counts()
    numbers=[int(largest.get('Microglia',0)),int(largest.get('Astrocyte',0)),int(largest.drop(['Microglia','Astrocyte'],errors='ignore').sum())]
    ax.barh(['Microglia','Astrocyte','All neurons combined'],numbers,color=[BLUE,ORANGE,GREY])
    for i,n in enumerate(numbers): ax.text(n+10,i,str(n),va='center')
    ax.set(xlim=(0,750),xlabel='Parents with this highest observed off-target',title='A  Related glia are the main specificity concern')
    ax.set_xticks([0,200,400,600])
    ax.invert_yaxis()
    ax=axs[0,1]
    depths=pseudo.groupby('analysis_group').retained_fragment_records.sum()/1e6
    names=['MG','Astrocyte','Microglia','Rod','Cone','OFF_cone_bipolar']
    ax.barh([n.replace('_',' ') for n in names],depths.loc[names],color=[BLUE,ORANGE,ORANGE,GREY,GREY,GREY])
    for i,n in enumerate(depths.loc[names]): ax.text(n+4,i,f'{n:.1f}',va='center',fontsize=9)
    ax.set(xlim=(0,445),xlabel='Retained pooled fragment records (millions)',title='B  Priority glial off-targets have limited depth')
    ax.invert_yaxis()
    ax=axs[1,0]
    names=['Repeat overlap ≥50%','k50 reference Umap mean <0.8','Parent orientation count gap >1','Attribution QC failure']
    values=[int(parents.repeat_overlap_fraction.ge(.5).sum()),int(parents.umap_k50_parent_mean.lt(.8).sum()),int(parents.parent_count_orientation_flag.sum()),int((~parents.attribution_qc_pass).sum())]
    ax.barh(names,values,color=[GREY,GREY,ORANGE,ORANGE])
    for i,n in enumerate(values): ax.text(n+5,i,str(n),va='center')
    ax.set(xlim=(0,320),xlabel='Flagged parents out of 1,000 (categories overlap)',title='C  Locus and orientation flags need review')
    ax.invert_yaxis()
    ax=axs[1,1];ax.axis('off')
    ax.text(0,1,'D  Missing evidence before a large order',weight='bold',fontsize=11,va='top')
    lines=[('Variant off-target accessibility','Not predicted: no off-target models'),
        ('Variant effects across training seeds','Not evaluated; only aggregate seed metrics'),
        ('Reporter activity / cell specificity','Not experimentally tested'),
        ('Unseen biological donors','Donor consistency only; no held-out-donor fit')]
    for i,(title,detail) in enumerate(lines):
        y=.82-i*.19
        ax.text(0,y,title,weight='bold',fontsize=10,color=RED,va='top')
        ax.text(0,y-.065,detail,fontsize=9,va='top')
    finish(fig,'04_open_questions','4. Evidence supports pilot planning, not a demonstrated MG-specific optimized library',
        'Before a large specificity-focused order: review representative triplets, quantify cross-seed edit stability, and assess glial off-target modeling.\nA limited pilot can test remaining uncertainty once assay/vector, controls, insert context and budget are specified.')
    metrics=dict(parents=len(parents),gain_variants=int(edits.sequence_role.eq('gain').sum()),loss_controls=int(edits.sequence_role.eq('loss_control').sum()),
        all_glial_lodo_fdr_pass=int(parents.maximum_glial_lodo_FDR.lt(.05).sum()),
        strongest_offtarget_microglia_or_astrocyte=numbers[0]+numbers[1],
        parent_flags={'repeat_ge_50pct':values[0],'k50_low':values[1],'orientation_gap_gt_1':values[2],'attribution_fail':values[3]},
        no_experimental_validation=True,no_offtarget_variant_models=True,
        decision='GO for pilot planning with explicit limitations; insufficient evidence for a large specificity-optimized synthesis order')
    (OUT/'evidence_summary.json').write_text(json.dumps(metrics,indent=2)+'\n')
    provenance={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(set(SOURCES))}
    (OUT/'source_checksums.json').write_text(json.dumps(provenance,indent=2)+'\n')
    titles=['Held-out model prediction and donor consistency','Natural-parent specificity and donor sensitivity','Variant screening and attribution consistency','Unresolved specificity and experimental risks']
    pngs=sorted(OUT.glob('0*.png'))
    intro='<h1>MG regulatory-element design: evidence before experiments</h1><p><b>Recommendation:</b> proceed with pilot planning; the current evidence does not establish that the optimized variants retain MG specificity or work in a reporter.</p><p>Existing results only. No new training or experimental validation was performed for this review. Generated for the October 14, 2026 project discussion.</p>'
    sections=''.join('<section><h2>'+title+'</h2><img alt="'+title+'" src="data:image/png;base64,'+base64.b64encode(path.read_bytes()).decode()+'"></section>' for title,path in zip(titles,pngs))
    html='<!doctype html><html lang="en"><meta charset="utf-8"><title>MG design evidence review</title><style>body{font:16px system-ui;max-width:1200px;margin:30px auto;padding:0 20px;color:#24313c}img{width:100%;height:auto}section{break-before:page;margin-top:40px}h1,h2{line-height:1.25}@media print{section{break-inside:avoid}body{margin:0}}</style>'+intro+sections+'<p>Source data, SHA256 manifests, numerical exports and the detailed interpretation accompany this report in the repository.</p></html>'
    (OUT/'evidence_review.html').write_text(html)
    print(json.dumps(metrics,indent=2))

if __name__=='__main__':
    main()
