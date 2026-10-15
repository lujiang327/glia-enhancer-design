#!/usr/bin/env python3
"""Verify the shareable track package and render predefined candidate loci."""
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
from prepare_browser_tracks import digest


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--tracks-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    for line in (a.tracks_dir/'checksums.sha256').read_text().splitlines():
        sha,name=line.split()
        if digest(a.tracks_dir/name)!=sha: raise ValueError('Checksum mismatch: '+name)
    manifest=json.loads((a.tracks_dir/'browser_track_manifest.json').read_text())
    candidates=pd.read_csv(a.tracks_dir/'candidate_index.tsv',sep='\t')
    bed=pd.read_csv(a.tracks_dir/'MG_top1000_parents.hg38.bed',sep='\t',header=None)
    assert len(bed)==len(candidates)==1000
    expected=candidates.sort_values(['chrom','start'])
    assert np.array_equal(bed.iloc[:,:3].to_numpy(),expected[['chrom','start','end']].to_numpy())
    sizes=dict((c,int(n)) for c,n in (line.split() for line in (a.tracks_dir/'hg38.canonical.chrom.sizes').read_text().splitlines()))
    handles=[]
    try:
        for track in manifest['tracks']:
            bw=pyBigWig.open(str(a.tracks_dir/track['file']))
            handles.append(bw)
            assert bw.chroms()==sizes
            assert np.isclose(bw.header()['sumData'],1e6,rtol=2e-5)
        ranks=[1,200,400,600,800,1000]
        colors=['black','#2878a5','#d88132','#659d62','#9c6cb0']
        fig,axes=plt.subplots(2,3,figsize=(13,7),sharey=True)
        for ax,rank in zip(axes.flat,ranks):
            row=candidates[candidates.candidate_rank.eq(rank)].iloc[0]
            left=max(0,int((row.start-1000)//100*100))
            right=min(sizes[row.chrom],int(np.ceil((row.end+1000)/100)*100))
            x=np.arange(left,right)-(row.start+250)
            for i,(bw,color,track) in enumerate(zip(handles,colors,manifest['tracks'])):
                values=np.asarray(bw.values(row.chrom,left,right),dtype=float)
                if not np.isfinite(values).all(): raise ValueError('Unexpected missing bin in preview')
                ax.step(x,values,where='mid',color=color,lw=1.4 if i==0 else .8,alpha=.85 if i==0 else .65,
                    label=track['file'].split('.IPM')[0].replace('MG_',''))
            ax.axvspan(-250,250,color='grey',alpha=.12)
            ax.set_title(f'Rank {rank}: {row.peak_id}\n{row.chrom}:{row.start+1}-{row.end}',fontsize=9)
            ax.set_xlim(x.min(),x.max())
            ax.spines[['top','right']].set_visible(False)
            ax.set_xlabel('Position relative to parent center (bp)',fontsize=9)
        for ax in axes[:,0]: ax.set_ylabel('Observed insertion density\n(IPM/base; 100bp mean)',fontsize=9)
        fig.suptitle('Observed MG ATAC signal: pooled and four donors',fontsize=15,weight='bold',y=.99)
        h,l=axes[0,0].get_legend_handles_labels()
        fig.legend(h,l,ncol=5,loc='upper center',bbox_to_anchor=(.5,.945),frameon=False)
        fig.text(.04,.02,'Predefined ranks 1, 200, 400, 600, 800, 1000; not selected for visual peak quality. Shading marks the natural 500bp parent.\nAll panels share the same y-scale. These are observed MG tracks only; off-target specificity and mutant activity cannot be inferred.',fontsize=9)
        fig.tight_layout(rect=[.02,.1,.99,.86],h_pad=2)
        a.output_dir.mkdir(parents=True,exist_ok=True)
        for ext in ['png','svg']: fig.savefig(a.output_dir/('MG_candidate_track_preview.'+ext),dpi=180)
        plt.close(fig)
    finally:
        for bw in handles: bw.close()
    (a.output_dir/'browser_track_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (a.output_dir/'browser_review_summary.json').write_text(json.dumps(dict(status='PASS_MG_BROWSER_TRACK_INTEGRITY',
        tracks=len(handles),candidates=len(candidates),preview_ranks=ranks,
        total_bigwig_bytes=sum(t['bytes'] for t in manifest['tracks']),
        off_target_tracks_included=False,all_bigwig_integrals_approximately_one_million=True),indent=2)+'\n')
    print((a.output_dir/'browser_review_summary.json').read_text())

if __name__=='__main__':main()
