#!/usr/bin/env python3
"""Render figures directly from measured QC tables; no imputed metrics."""
import csv
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR',str(Path('data/intermediate/matplotlib').resolve()))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr


def rows(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))


def run():
    out=Path('reports/qc');figdir=Path('reports/figures');figdir.mkdir(exist_ok=True)
    donors=rows(out/'donor_summary.tsv');names=[r['donor'] for r in donors]
    bins=rows(out/'donor_50kb_bins.tsv');counts=np.array([[float(r[d]) for d in names] for r in bins]);counts=counts[counts.sum(1)>0]
    cpm=counts/counts.sum(0)*1e6;logcpm=np.log1p(cpm);corr=np.corrcoef(logcpm.T);rho=spearmanr(counts,axis=0).statistic
    with (out/'donor_correlations.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['donor_a','donor_b','pearson_log1p_CPM_50kb','spearman_50kb'])
        for i,a in enumerate(names):
            for j,b in enumerate(names):w.writerow([a,b,corr[i,j],rho[i,j]])
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
    colors=['#287271','#457b9d','#ad6d42','#74618c']
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    ax=axs[0,0];v=[int(r['unambiguous_MG_cells']) for r in donors];bars=ax.bar(names,v,color=colors);ax.bar_label(bars);ax.set_ylim(0,max(v)*1.15);ax.set_title('A  Retained Müller glia cohort');ax.set_ylabel('Unambiguously mapped cells')
    ax=axs[0,1];v=[int(r['retained_fragment_records'])/1e6 for r in donors];bars=ax.bar(names,v,color=colors);ax.bar_label(bars,fmt='%.2f');ax.set_ylim(0,max(v)*1.15);ax.set_title('B  Donor pseudobulk depth');ax.set_ylabel('Filtered fragment records (millions)')
    ts=rows(out/'donor_tss_profiles.tsv');x=np.array([int(r['offset']) for r in ts]);ax=axs[1,0]
    for d,color in zip(names,colors):
        y=np.array([float(r[d]) for r in ts]);base=np.r_[y[:100],y[-100:]].mean();sm=np.convolve(y/base,np.ones(25)/25,'valid');ax.plot(x[12:-12],sm,label=d,color=color)
    ax.set_title('C  Pooled donor TSS profiles');ax.set_xlabel('Distance to RefGene TSS (bp)');ax.set_ylabel('Signal / flank mean (25-bp smooth)');ax.legend(fontsize=8)
    ax=axs[1,1];im=ax.imshow(corr,vmin=0,vmax=1,cmap='YlGnBu');ax.set_xticks(range(4),names);ax.set_yticks(range(4),names);ax.set_title('D  Donor consistency: 50-kb QC bins')
    for i in range(4):
        for j in range(4):ax.text(j,i,f'{corr[i,j]:.3f}',ha='center',va='center',color='white' if corr[i,j]>.65 else 'black')
    fig.colorbar(im,ax=ax,label='Pearson r, log1p CPM');fig.suptitle('GSE196235 • donor-resolved subset; both eyes pooled per donor',fontsize=14)
    fig.savefig(figdir/'donor_qc.png',dpi=160);fig.savefig(figdir/'donor_qc.svg');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(13,5),layout='constrained');ax=axs[0]
    lengths=rows(out/'fragment_lengths.tsv');x=np.array([int(r['fragment_lengths']) for r in lengths]);y=np.array([int(r['records']) for r in lengths]);keep=x<=1000;ax.plot(x[keep],y[keep]/y.sum()*100,color='#287271');ax.set_title('A  Full-pool fragment length distribution');ax.set_xlabel('Fragment length (bp)');ax.set_ylabel('% of fragment records per bp')
    markers=rows(out/'marker_summary.tsv');lookup={(r['group'],r['gene']):r for r in markers};genes=['GLUL','RLBP1','SLC1A3','SOX9','AQP4','GFAP','PDE6A','RHO','ARR3','GRM6','GAD1','C1QA'];groups=['MG_'+d for d in names]+['non_MG_barcode_background'];ax=axs[1]
    for j,g in enumerate(groups):
        v=[lookup[(g,m)] for m in genes];sc=ax.scatter(range(len(genes)),[j]*len(genes),s=[5+130*float(r['fraction_detected']) for r in v],c=[float(r['mean_log1p_CP10k']) for r in v],vmin=0,vmax=5,cmap='viridis')
    ax.set_xticks(range(len(genes)),genes,rotation=65,ha='right');ax.set_yticks(range(len(groups)),names+['Other barcode background']);ax.invert_yaxis();ax.set_title('B  Published-label RNA marker check');fig.colorbar(sc,ax=ax,label='Mean log1p(CP10k)');ax.text(0,-.22,'Dot area ∝ detected fraction; background is not a pure cell type',transform=ax.transAxes,fontsize=8)
    fig.suptitle('GSE196235 • marker support with donor-dependent neuronal RNA',fontsize=14);fig.savefig(figdir/'markers_and_lengths.png',dpi=160);fig.savefig(figdir/'markers_and_lengths.svg');plt.close(fig)
    print('Saved four figure files and measured donor correlations.')


if __name__=='__main__':run()
