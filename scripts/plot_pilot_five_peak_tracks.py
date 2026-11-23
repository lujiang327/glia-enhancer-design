#!/usr/bin/env python3
"""Plot the five user-specified pilot parents from verified observed BigWigs."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import pyBigWig

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reports/browser_tracks/pilot_five_candidates'
PDF=ROOT/'output/pdf/MG_five_pilot_candidate_peak_tracks.pdf'
TRACKS=ROOT/'data/intermediate/browser_tracks/retina_hg38_review_100bp'
OUT.mkdir(parents=True,exist_ok=True);PDF.parent.mkdir(parents=True,exist_ok=True)
RANKS=[1,47,140,25,779]
GROUPS=['MG','Astrocyte','Microglia','Rod','Cone','ON_cone_bipolar','OFF_cone_bipolar',
        'Rod_bipolar','GABA_amacrine','Gly_amacrine','AII_amacrine','Horizontal','RGC']
DONORS=['LGS1','LGS2','LGS3','LVG1']
LABELS={'MG':'Müller glia','Astrocyte':'Astrocytes','Microglia':'Microglia',
    'Rod':'Rods','Cone':'Cones','ON_cone_bipolar':'ON cone bipolar',
    'OFF_cone_bipolar':'OFF cone bipolar','Rod_bipolar':'Rod bipolar',
    'GABA_amacrine':'GABA amacrine','Gly_amacrine':'Gly amacrine',
    'AII_amacrine':'AII amacrine','Horizontal':'Horizontal','RGC':'Ganglion cells'}
COLORS={'MG':'#192B26','Astrocyte':'#D47D31','Microglia':'#287DB2'}


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


source=ROOT/'reports/phase3/candidate_ranking/top_1000_candidates.tsv.gz'
parents=pd.read_csv(source,sep='\t').set_index('candidate_rank').loc[RANKS].reset_index()
audit=pd.read_csv(ROOT/'reports/browser_tracks/parent_specificity_audit/browser_signal.tsv',sep='\t').set_index('candidate_rank')
donor_audit=pd.read_csv(ROOT/'reports/browser_tracks/parent_specificity_audit/donor_browser_signal.tsv',sep='\t').set_index('candidate_rank')
manifest=pd.read_csv(TRACKS/'track_index.tsv',sep='\t').set_index('file')
assert parents.candidate_rank.tolist()==RANKS and parents.end.sub(parents.start).eq(500).all()
names=[g+'_pooled' for g in GROUPS]+[g+'_'+d for d in DONORS for g in GROUPS[:3]]
values=np.zeros((5,len(names),2500))
hashes={}
for j,name in enumerate(names):
    p=TRACKS/(name+'.IPM_per_base.100bp.bw')
    hashes[p.name]=sha(p)
    assert hashes[p.name]==manifest.loc[p.name,'sha256']
    with pyBigWig.open(str(p)) as b:
        for i,r in parents.iterrows():
            left,right=int(r.start)-1000,int(r.end)+1000
            assert left>=0 and right<=b.chroms(r.chrom)
            v=np.asarray(b.values(r.chrom,left,right),dtype=float)
            assert len(v)==2500 and np.isfinite(v).all() and (v>=0).all()
            values[i,j]=v
means=values[:,:,1000:1500].mean(axis=2)
for i,r in parents.iterrows():
    for j,g in enumerate(GROUPS):
        assert np.isclose(means[i,j],audit.loc[r.candidate_rank,g],rtol=1e-12,atol=1e-15)
    for j,n in enumerate(names[13:],13):
        assert np.isclose(means[i,j],donor_audit.loc[r.candidate_rank,n],rtol=1e-12,atol=1e-15)
parents['IGV_parent_locus']=parents.chrom+':'+(parents.start+1).astype(str)+'-'+parents.end.astype(str)
parents['browser_MG_mean_IPM_per_base']=means[:,0]
parents['browser_MG_to_maxoff_ratio']=means[:,0]/means[:,1:13].max(axis=1)
parents['browser_strongest_offtarget']=[GROUPS[i+1] for i in means[:,1:13].argmax(axis=1)]
parents['browser_donors_MG_above_both_glia']=np.sum(np.stack([means[:,13+k*3]>means[:,14+k*3:16+k*3].max(axis=1) for k in range(4)]),axis=0)
cols=['candidate_rank','peak_id','IGV_parent_locus','genomic_class','nearest_gene','muller_mean_cpm',
      'browser_MG_mean_IPM_per_base','browser_MG_to_maxoff_ratio','browser_strongest_offtarget',
      'primary_logFC','primary_FDR','primary_donors_positive','browser_donors_MG_above_both_glia']
parents[cols].to_csv(OUT/'five_candidate_metrics.tsv',sep='\t',index=False)
signal=parents[['candidate_rank','peak_id','IGV_parent_locus']].copy()
for j,n in enumerate(names):signal[n+'_mean_IPM_per_base']=means[:,j]
signal.to_csv(OUT/'five_candidate_track_means.tsv',sep='\t',index=False)
np.savez_compressed(OUT/'displayed_track_values.npz',values=values,
    relative_positions=np.arange(2500)-1250,track_names=np.asarray(names),candidate_ranks=np.asarray(RANKS))
with (OUT/'five_candidates.hg38.bed').open('w') as f:
    for _,r in parents.iterrows():f.write(f'{r.chrom}\t{r.start}\t{r.end}\tP{r.candidate_rank}_{r.nearest_gene}\n')

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
x=np.arange(2500)-1250
pooled_max=float(values[:,:13].max())*1.08
donor_max=float(values[:,13:].max())*1.08


def draw(page,zoom=False,donors=False):
    track_indices=list(range(13,len(names))) if donors else list(range(13))
    ymax=donor_max if donors else (float(values[:,:13,1000:1500].max())*1.08 if zoom else pooled_max)
    fig=plt.figure(figsize=(19.5,12.6))
    if donors:
        title='Do the MG-enriched peaks recur across donors?'
        subtitle='Four donors; MG, astrocytes and microglia for each donor. One y-scale across all 60 donor panels.'
    elif zoom:
        title='The exact 500-bp sequences proposed for cloning'
        subtitle='Zoom of the same observed tracks. One shared y-scale across all five candidates and all 13 cell types.'
    else:
        title='Five pilot candidates: observed retinal ATAC peak tracks'
        subtitle='Same 2.5-kb window and one shared y-scale across all five candidates and all 13 cell types.'
    fig.text(.035,.962,title,fontsize=20,weight='bold')
    fig.text(.035,.928,subtitle,fontsize=11)
    fig.text(.035,.903,'MG: dark green   |   Astrocytes: orange   |   Microglia: blue   |   Other retinal cells: gray   |   Shading: exact 500-bp parent',fontsize=10,color='#425466')
    gs=fig.add_gridspec(len(track_indices),5,left=.125,right=.98,bottom=.14,top=.77,hspace=.16,wspace=.19)
    if not zoom:
        header_y=.866
    else:header_y=.866
    for c,r in parents.iterrows():
        gene=str(r.nearest_gene);kind={'distal_intergenic':'distal','promoter':'promoter-associated','intronic':'intronic'}.get(r.genomic_class,r.genomic_class)
        hdr=(f'P{int(r.candidate_rank)}  |  {gene}  ({kind})\n{r.IGV_parent_locus}\n'
             f'MG mean {r.browser_MG_mean_IPM_per_base:.4f}  |  MG/max off-target {r.browser_MG_to_maxoff_ratio:.2f}x\n'
             f'MG CPM {r.muller_mean_cpm:.1f}  |  donor support {int(r.browser_donors_MG_above_both_glia)}/4')
        pos=.125+c*((.98-.125)/5)+(.98-.125)/10
        fig.text(pos,header_y,hdr,ha='center',va='top',fontsize=8.8,linespacing=1.4)
        for ri,j in enumerate(track_indices):
            name=names[j];g=name.split('_')[0] if donors else GROUPS[j]
            ax=fig.add_subplot(gs[ri,c]);v=values[c,j]
            ax.axvspan(-250,250,color='#DCE5ED',alpha=.65,zorder=0)
            # These are already 100bp-binned values; no smoothing or autoscaling.
            ax.fill_between(x,v,step='mid',color=COLORS.get(g,'#8D99A6'),linewidth=0,zorder=2)
            ax.set_xlim(-250 if zoom else -1250,250 if zoom else 1250);ax.set_ylim(0,ymax)
            ax.spines[['left','top','right']].set_visible(False)
            ax.spines['bottom'].set_color('#DCE2E7');ax.set_yticks([])
            if c==0:
                label=(name.split('_',1)[1]+'  '+LABELS[g]) if donors else LABELS[g]
                ax.set_ylabel(label,rotation=0,ha='right',va='center',fontsize=8.5,labelpad=9,
                    color=COLORS.get(g,'#425466'),weight='bold' if g=='MG' else 'normal')
            if donors and ri%3==0 and ri>0:ax.axhline(ymax,color='#BBC6D0',lw=.7)
            if ri==0:
                ax.spines['left'].set_visible(True);ax.spines['left'].set_color('#8795A0')
                ax.set_yticks([0,ymax]);ax.set_yticklabels(['0',f'{ymax:.3f}'],fontsize=7);ax.tick_params(axis='y',length=2,pad=2)
            if ri==len(track_indices)-1:
                ax.set_xticks([-250,0,250] if zoom else [-1000,0,1000]);ax.tick_params(axis='x',labelsize=8)
                ax.set_xlabel('bp from parent center',fontsize=9)
            else:ax.set_xticks([])
    fig.text(.035,.092,f'Every panel on this page uses y = 0 to {ymax:.5f} IPM/base. Original 100-bp binning; no smoothing, peak-height matching or individual autoscaling.',fontsize=10)
    fig.text(.035,.065,'Compare the dark-green signal inside the shaded insert, then residual off-target signal. Peaks outside the shading are not part of the cloned fragment.',fontsize=10)
    fig.text(.035,.040,'Track density: observed fragment-end insertions per million insertions per base. Header CPM: equal-donor TMM accessibility; these normalizations are different.',fontsize=9,color='#425466')
    fig.text(.035,.020,f'Natural-locus accessibility is not reporter activity or experimentally proven MG specificity. Nearest genes are labels, not verified targets.   |   {page}/3',fontsize=9,color='#425466')
    label='03_donor_glial_tracks' if donors else ('02_parent_500bp_zoom' if zoom else '01_pooled_retinal_tracks')
    fig.savefig(OUT/(label+'.png'),dpi=175,facecolor='white')
    return fig,ymax


with PdfPages(PDF) as pdf:
    pdf.infodict().update(Title='Five MG pilot candidates: observed peak tracks',Author='Glia regulatory-element project')
    ymaxes=[]
    for page,kwargs in [(1,{}),(2,dict(zoom=True)),(3,dict(donors=True))]:
        fig,ymax=draw(page,**kwargs);pdf.savefig(fig,facecolor='white');plt.close(fig);ymaxes.append(ymax)
summary=dict(date='2026-11-16',candidate_ranks=RANKS,cell_types=GROUPS,donors=DONORS,
    pages=['2.5kb pooled all retinal cells','exact 500bp pooled zoom','2.5kb donor glial comparisons'],
    global_ymax_per_page=ymaxes,normalization='Original global IPM/base, 100bp bins; no rescaling',
    source_tracks_sha256=hashes,source_candidates_sha256=sha(source),
    all_source_track_hashes_match_manifest=True,all_500bp_means_match_previous_audit=True,
    candidate_selection='Five user-specified candidates; no track-shape-based re-selection',
    limitations=['Already 100bp binned, not base-resolution insertion signal','Observed natural-locus accessibility, not AAV reporter activity',
                 'Pooled and donor pages have separately stated scales; scales identical across panels within a page'],
    pdf=str(PDF.relative_to(ROOT)))
(OUT/'provenance_and_validation.json').write_text(json.dumps(summary,indent=2)+'\n')
(OUT/'README.md').write_text('''# Five pilot candidate peak tracks

P1, P47, P140, P25 and P779 are shown in that order. These are the original
candidate ranks, not a new ranking.

Open `output/pdf/MG_five_pilot_candidate_peak_tracks.pdf` from the repository root.
The same pages are saved as PNG files in this folder.

1. **Pooled retinal tracks:** 2.5 kb around each parent, all13 cell types. Compare
   strength across columns and off-target signal down each column; all panels
   share one y-scale.
2. **Exact500-bp zoom:** focuses on the sequences to be cloned. This page has its
   own shared y-scale, explicitly stated; do not compare its heights directly to
   page1 without accounting for the scale.
3. **Donor glial tracks:** each donor's MG, astrocyte and microglia profiles. All
   panels share one donor-page scale. Compare within each donor and candidate.

Shading marks the cloned500 bp. Signal outside that interval is genomic context,
not included in the insert. Tracks retain the original global IPM/base scaling
and100-bp bins, without smoothing or individual autoscaling. Header CPM comes
from equal-donor TMM accessibility, a different normalization. Neither measure
establishes reporter activity, delivery or experimental cell-type specificity.

`five_candidate_metrics.tsv` records metrics; `five_candidate_track_means.tsv`
records pooled/donor means; `displayed_track_values.npz` preserves exact plotted
samples. The BED file provides the five hg38 intervals. Source hashes, settings
and cross-checks against the earlier audit are in `provenance_and_validation.json`.

Reproduce with `MPLCONFIGDIR=/tmp/glia-matplotlib .venv/bin/python scripts/plot_pilot_five_peak_tracks.py`.
''')
print(json.dumps(dict(pdf=str(PDF),ranks=RANKS,pages=3,ymaxes=ymaxes,source_tracks_verified=len(hashes)),indent=2))
