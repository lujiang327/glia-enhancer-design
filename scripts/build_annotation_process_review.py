#!/usr/bin/env python3
"""Document inherited retinal annotation and our measured processing/QC evidence."""
import hashlib
import gzip
import json
from collections import Counter
from pathlib import Path
import shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reports/annotation_process_review'
PDF=ROOT/'output/pdf/MG_multiome_processing_annotation_review.pdf'
DONORS=['LGS1','LGS2','LGS3','LVG1']
GROUPS=['MG','Astrocyte','Microglia','Rod','Cone','ON_cone_bipolar','OFF_cone_bipolar',
        'Rod_bipolar','GABA_amacrine','Gly_amacrine','AII_amacrine','Horizontal','RGC']
SOURCES={}


def source(name):
    p=ROOT/name
    SOURCES[name]=hashlib.sha256(p.read_bytes()).hexdigest()
    return p


def table(name):return pd.read_csv(source(name),sep='\t')


def plots():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    # Distinguish the published analysis from our reuse; no invented clustering.
    fig,ax=plt.subplots(figsize=(10,7.7));ax.axis('off');ax.set_xlim(0,1);ax.set_ylim(0,1)
    lanes=[(.025,'Original study: reported workflow','#F6EBD9',[
        ('1. Eight retinas, four donors','Paired 10x multiome RNA + ATAC\nCell Ranger ARC 2.0.0; hg38'),
        ('2. Joint cell QC and doublets','RNA UMI >200 and <50,000; mito <1%; ribo <5%\nATAC TSS >6; fragments >2,500; ArchR doublets'),
        ('3. RNA clustering','Seurat normalization; 2,000 variable genes; PCA\nHarmony by sample; first 20 dimensions; resolution 0.5'),
        ('4. Marker-based annotation','Remove mixed-marker / marker-poor clusters; recluster\n22 clusters grouped into 13 cell types'),
        ('5. Transfer paired labels to ATAC','RNA and ATAC linked by nucleus barcode\nPublish RNA matrix and cell-type fragment pools')]),
        (.535,'Our project: completed reuse and QC','#E5EFF5',[
        ('1. Reuse processed GSE196235','51,645 RNA columns; 13 labelled ATAC pools\nNo new alignment, RNA clustering or manual relabeling'),
        ('2. Recover donor identities','Use archived sample-prefix mapping\nKeep unique RNA-barcode matches; exclude ambiguity'),
        ('3. Check MG annotation and ATAC quality','MG RNA markers by donor; TSS proxy / FRiP\nFragment depth, length and donor consistency'),
        ('4. Preserve four biological donors','3,434 MG cells; 85,979,585 retained MG fragments\n52 donor-by-cell-type ATAC pseudobulks'),
        ('5. Model, then compare accessibility','Genome-wide ChromBPNet training and QC first\nDonor-blocked DA; candidates; observed browser tracks')])]
    for x,title,fill,steps in lanes:
        ax.text(x+.22,.975,title,ha='center',va='top',fontsize=12,weight='bold')
        for i,(head,body) in enumerate(steps):
            y=.77-i*.18
            ax.add_patch(FancyBboxPatch((x,y),.44,.145,boxstyle='round,pad=0.009',facecolor=fill,edgecolor='#B7C7D2'))
            ax.text(x+.015,y+.115,head,fontsize=10,weight='bold',va='top')
            ax.text(x+.015,y+.073,body,fontsize=8.5,va='top',linespacing=1.5)
            if i<4:ax.annotate('',xy=(x+.22,y-.025),xytext=(x+.22,y-.006),arrowprops=dict(arrowstyle='->',color='#637787'))
    fig.tight_layout(pad=.1);fig.savefig(OUT/'01_processing_workflow.png',dpi=190);plt.close(fig)
    markers=table('reports/qc/marker_summary.tsv')
    genes=['GLUL','RLBP1','SLC1A3','SOX9','AQP4','GFAP','PDE6A','RHO','ARR3','GRM6','GAD1','C1QA']
    rows=['MG_all']+['MG_'+d for d in DONORS]+['non_MG_barcode_background']
    fig,ax=plt.subplots(figsize=(10,4.8))
    for i,g in enumerate(rows):
        for j,gene in enumerate(genes):
            r=markers[markers.group.eq(g)&markers.gene.eq(gene)].iloc[0]
            im=ax.scatter(j,i,s=240*r.fraction_detected,c=[r.mean_log1p_CP10k],cmap='viridis',vmin=0,vmax=3.5,edgecolors='#808080',linewidths=.3)
    labels=[f'All MG (n=3,434)']+[f'{d} MG (n={int(markers.loc[markers.group.eq("MG_"+d),"n_cells"].iloc[0]):,})' for d in DONORS]+['Other RNA barcodes (mixed types)']
    ax.set(xticks=range(12),xticklabels=genes,yticks=range(6),yticklabels=labels,xlim=(-.7,11.7),ylim=(5.6,-.7))
    ax.tick_params(axis='x',rotation=50)
    ax.axvline(5.5,color='#ADB9C2',ls='--',lw=.8)
    ax.set_title('Existing RNA marker check of published MG-labelled nuclei',weight='bold',pad=20)
    cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.025);cb.set_label('Mean log1p(CP10k), no per-gene scaling')
    handles=[ax.scatter([],[],s=240*v,color='#999999') for v in [.25,.5,1]]
    ax.legend(handles,['25%','50%','100%'],title='Cells detected',ncol=3,frameon=False,loc='upper center',bbox_to_anchor=(.5,-.33))
    fig.tight_layout(rect=[0,.04,1,.96]);fig.savefig(OUT/'03_MG_RNA_marker_review.png',dpi=190);plt.close(fig)
    markers.to_csv(OUT/'RNA_marker_evidence.tsv',sep='\t',index=False)
    cells=table('reports/phase3/donor_celltype_audit.tsv')
    counts=cells.pivot(index='analysis_group',columns='donor',values='unambiguous_cells').loc[GROUPS,DONORS]
    counts.to_csv(OUT/'donor_celltype_counts.tsv',sep='\t')
    pca=table('reports/phase3/count_matrix_qc/pca_scores.tsv')
    q=json.loads(source('reports/phase3/count_matrix_qc/count_matrix_qc.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(10,7.2),gridspec_kw={'width_ratios':[.9,1.2]})
    ax=axes[0];im=ax.imshow(np.log10(counts.to_numpy()),cmap='Blues',aspect='auto',vmin=1,vmax=4.2)
    for i in range(13):
        for j in range(4):
            n=counts.iloc[i,j]
            ax.text(j,i,f'{n:,}',ha='center',va='center',fontsize=8,color='white' if n>2000 else '#24313C',weight='bold' if n<50 else 'normal')
            if n<50:ax.add_patch(plt.Rectangle((j-.5,i-.5),1,1,fill=False,edgecolor='#C85E26',lw=2))
    ax.set(xticks=range(4),xticklabels=DONORS,yticks=range(13),yticklabels=[g.replace('_',' ') for g in GROUPS],title='A  Donor-resolved cells after barcode exclusions')
    ax.tick_params(labelsize=8)
    ax=axes[1];palette=plt.get_cmap('tab20');shapes=['o','s','^','D']
    for i,g in enumerate(GROUPS):
        for j,d in enumerate(DONORS):
            r=pca[pca.analysis_group.eq(g)&pca.donor.eq(d)].iloc[0]
            ax.scatter(r.PC1,r.PC2,color=palette(i),marker=shapes[j],s=50,label=g.replace('_',' ') if j==0 else None,edgecolors='white',linewidths=.3)
    ax.set(title='B  ATAC pseudobulk PCA (52 samples)',xlabel=f'PC1 ({q["pca_associations"][0]["explained_variance_fraction"]:.1%})',ylabel=f'PC2 ({q["pca_associations"][1]["explained_variance_fraction"]:.1%})')
    ax.legend(loc='upper center',bbox_to_anchor=(.50,-.16),ncol=2,fontsize=7,frameon=False)
    fig.text(.04,.025,'Orange outlines: fewer than 50 cells. PCA is a sample-level accessibility diagnostic, not a new single-cell UMAP or independent label assignment.',fontsize=9)
    fig.tight_layout(rect=[0,.15,1,.96],w_pad=2.5);fig.savefig(OUT/'04_donor_celltype_ATAC_review.png',dpi=190);plt.close(fig)
    return markers,counts,q


def build_pdf(markers,counts,q):
    styles=getSampleStyleSheet()
    styles.add(ParagraphStyle(name='ATitle',fontName='Helvetica-Bold',fontSize=21,leading=25,spaceAfter=9,textColor=colors.HexColor('#243746')))
    styles.add(ParagraphStyle(name='ASub',fontName='Helvetica-Bold',fontSize=13,leading=17,spaceBefore=8,spaceAfter=5,textColor=colors.HexColor('#216B8E')))
    styles.add(ParagraphStyle(name='ABody',fontName='Helvetica',fontSize=10.4,leading=14.4,spaceAfter=7))
    styles.add(ParagraphStyle(name='ASmall',fontName='Helvetica',fontSize=8.8,leading=11.5,spaceAfter=5,textColor=colors.HexColor('#566775')))
    def p(s,style='ABody'):return Paragraph(s,styles[style])
    def img(name,width=516):
        f=OUT/name;w,h=PILImage.open(f).size
        return Image(str(f),width=width,height=width*h/w)
    def tab(rows,widths):
        t=Table([[p(str(x),'ASmall') for x in r] for r in rows],colWidths=widths,hAlign='LEFT')
        t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#E5EFF5')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#AFC4D0')),('BOTTOMPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),6)]))
        return t
    def footer(c,d):
        c.setFont('Helvetica',8);c.setFillColor(colors.HexColor('#566775'))
        c.drawString(48,26,'Retinal MG project | Processing and annotation evidence | October 20, 2026')
        c.drawRightString(564,26,str(d.page))
    story=[p('How the retinal cells were processed and annotated','ATitle'),
        p('<b>Scope:</b> GSE196235 adult human retina. We retained published labels and independently checked aspects of the processed data. We have not rerun raw-read processing or reproduced the original single-cell clustering.'),
        img('01_processing_workflow.png'),
        p('The original study profiled 51,645 nuclei. Our donor-resolved analysis retains 48,422 cell records across 13 published types, including 3,434 MG nuclei. The two eyes of a donor are combined; they are not independent biological donors.','ASmall'),
        p('Sources: Wang et al., Cell Genomics 2022, doi:10.1016/j.xgen.2022.100164; archived RNA/ATAC analysis code, doi:10.5281/zenodo.6795162. Exact local source references accompany this packet.','ASmall'),
        PageBreak(),p('Original study: clusters and marker annotation','ATitle'),
        p('<b>This is the authors\' published figure, not a UMAP generated by our pipeline.</b> It shows their 22 RNA clusters grouped into 13 cell types and the marker-expression panel used to support those identities.'),
        img('02_published_Wang2022_Figure1.jpg',500),
        p('<b>Read B:</b> the published RNA UMAP. <b>Read D:</b> marker expression by published type. GLUL/RLBP1 support MG; PAX2/GFAP and C1QA/HLA-DRA help distinguish astrocytes and microglia. Marker expression can overlap between types.'),
        p('Wang et al. (2022), Figure 1, reproduced unchanged from PMC9584034. Copyright 2022 The Author(s), CC BY 4.0. https://pmc.ncbi.nlm.nih.gov/articles/PMC9584034/#fig1. This figure is contextual published evidence, not independent validation of our barcode mapping.','ASmall'),
        PageBreak(),p('Our RNA check: MG support and a contamination caveat','ATitle'),
        p('These plots use our existing matched-RNA marker summaries for the donor-resolved MG subset. Dot size is the fraction with at least one detected transcript; color is mean log1p expression after scaling each nucleus to 10,000 RNA counts.'),
        img('03_MG_RNA_marker_review.png'),
        p('Marker detection in the retained MG cohort','ASub')]
    rows=[['Group','GLUL','RLBP1','PDE6A','RHO']]
    for group in ['MG_all']+['MG_'+d for d in DONORS]:
        subset=markers[markers.group.eq(group)].set_index('gene')
        rows.append([group.replace('MG_all','All MG').replace('MG_','')]+[f'{subset.loc[g,"fraction_detected"]:.1%}' for g in ['GLUL','RLBP1','PDE6A','RHO']])
    story.extend([tab(rows,[136,95,95,95,95]),Spacer(1,8),
        p('<b>Support:</b> GLUL is detected in 97.1% and RLBP1 in 86.7% of the retained MG-labelled nuclei. Detection varies across donors.'),
        p('<b>Unresolved:</b> rod-associated PDE6A is detected in 55.5% and RHO in 63.7%. Ambient RNA or residual mixed nuclei are plausible explanations, but neither has been demonstrated or corrected. This check does not independently validate all 13 cell-type labels, and the mixed background is not a dedicated astrocyte/microglia comparison.'),
        PageBreak(),p('Our ATAC checks: donor coverage and sample structure','ATitle'),
        img('04_donor_celltype_ATAC_review.png'),
        p('<b>Coverage:</b> all 13 types are represented in all four donors. Orange boxes flag astrocyte LGS1 (12 cells), astrocyte LGS2 (46), microglia LVG1 (27), and RGC LGS1 (42).'),
        p('<b>ATAC structure:</b> the 52 pseudobulks separate by published cell type more strongly than by donor in the leading PCs. This supports broad coherence but cannot prove every nucleus is correctly annotated. Some related-glial samples correlate more strongly with same-donor MG than their own group.'),
        p('<b>Depth and quality:</b> 775,669,056 retained fragments across 52 pseudobulks; per-sample consensus FRiP 0.299-0.633. MG has 85,979,585 retained fragments. Its donor TSS proxies are 20.20-22.76, and diagnostic FRiP is 61.5%-64.5%. The pooled TSS proxy is not the authors\' per-cell ArchR score; FRiP changes with the peak set.','ASmall'),
        PageBreak(),p('What is verified, and what remains open','ATitle'),
        tab([
          ['Item','Evidence available','Scope / remaining limit'],
          ['Published annotation','Paper, original UMAP / dot plot, archived RNA and ATAC scripts','Inherited labels; no independent reproduction of the final 22-cluster object.'],
          ['Barcode to donor','Explicit eight-sample prefix map; ambiguous bare barcodes excluded','MG: 231 ambiguous bare barcodes excluded. Full donor counts for excluded cells remain unresolved.'],
          ['MG marker support','Matched RNA summaries by donor; GLUL/RLBP1 support identity','Rod RNA signal remains. No ambient correction or fresh doublet inference was performed.'],
          ['ATAC processing','Canonical hg38, valid coordinates, full-fragment blacklist filtering; donor-specific pseudobulks','Alignment / deduplication and original per-cell filtering inherited from processed release.'],
          ['Broader cell types','13 file labels, donor coverage, pseudobulk PCA and correlations','A fresh all-13-type RNA marker audit and global cross-pool cell-ID conflict check are not supplied here.']
        ],[105,205,206]),
        p('Two code-level reproducibility items','ASub'),
        p('The archived RNA script explicitly uses Harmony for initial clustering, but its later reclustering block omits the reduction argument. The archived ATAC snippet assigns RNA identities as positional vectors without a visible named-cell reorder assertion. These are items to verify with the final objects and execution logs, not proof that the released annotations are wrong.'),
        p('Recommended next check','ASub'),
        p('Obtain the full per-cell sample-prefixed identity table and final RNA object / embedding. Verify one-to-one RNA-to-ATAC cell joins, donor/eye mapping, and absence of conflicting labels across pools. Plot canonical markers separately for all 13 types, especially MG versus astrocytes/microglia; inspect rod-high MG cells and donor dependence before deciding whether a sensitivity subset is needed. Preserve published labels unless this reveals a clear error.'),
        p('<b>Current conclusion:</b> published annotation plus our QC supports continued exploratory analysis. We cannot yet certify that all upstream processing and annotation are issue-free. Strong model or DA results alone do not resolve that uncertainty.'),
        p('Primary sources: https://doi.org/10.1016/j.xgen.2022.100164; https://doi.org/10.5281/zenodo.6795162; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196235. Local scripts, tables, source hashes and detailed line references are in reports/annotation_process_review/README.md.','ASmall')])
    doc=SimpleDocTemplate(str(PDF),pagesize=letter,leftMargin=48,rightMargin=48,topMargin=38,bottomMargin=46,
                          title='Retinal multiome processing and cell-type annotation review',author='Glia regulatory-element project')
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    assert len(PdfReader(str(PDF)).pages)==5


def main():
    OUT.mkdir(parents=True,exist_ok=True);PDF.parent.mkdir(parents=True,exist_ok=True)
    image=source('data/intermediate/annotation_review_sources/Wang2022_Figure1.jpg')
    shutil.copy2(image,OUT/'02_published_Wang2022_Figure1.jpg')
    for name in ['data/metadata/author_code/seankwang-retina-multiome-651e181/scRNA',
                 'data/metadata/author_code/seankwang-retina-multiome-651e181/scATAC',
                 'config/samples.tsv','config/phase3_cell_types.tsv','reports/qc/fragment_audit.json',
                 'data/metadata/GSE196235_asset_inventory.json','reports/qc/marker_qc.json',
                 'reports/qc/donor_metadata.tsv','data/raw/GSM5866081_barcodes.tsv.gz',
                 'reports/qc/donor_summary.tsv','reports/qc/frip.tsv','reports/phase3/pseudobulk_qc.tsv',
                 'scripts/marker_qc.py','scripts/audit_fragments.py','scripts/audit_phase3_celltype_fragments.py']:
        source(name)
    markers,counts,q=plots()
    assert int(counts.to_numpy().sum())==48422 and int(counts.loc['MG'].sum())==3434
    ids=table('reports/qc/cells.unambiguous.tsv')
    assert len(ids)==3434 and ids.cell_id.is_unique and ids.bare_barcode.is_unique
    samples=table('config/samples.tsv').set_index('rna_prefix')
    for r in ids.itertuples():
        assert r.cell_id==f'{r.rna_prefix}_{r.bare_barcode}'
        assert r.donor==samples.loc[r.rna_prefix,'donor'] and r.eye==samples.loc[r.rna_prefix,'eye']
    with gzip.open(ROOT/'data/raw/GSM5866081_barcodes.tsv.gz','rt') as handle:
        rna=[line.strip() for line in handle if line.strip()]
    assert len(rna)==len(set(rna))==51645 and set(ids.cell_id)<=set(rna)
    prefix_counts=Counter(cell.split('_',1)[0] for cell in rna)
    assert set(prefix_counts)==set('12345678')
    sample_counts=table('config/samples.tsv')
    sample_counts['RNA_matrix_columns']=sample_counts.rna_prefix.astype(str).map(prefix_counts)
    sample_counts['retained_MG_nuclei']=sample_counts.rna_prefix.map(ids.rna_prefix.value_counts()).fillna(0).astype(int)
    assert sample_counts.retained_MG_nuclei.sum()==3434
    sample_counts.to_csv(OUT/'RNA_and_MG_sample_counts.tsv',sep='\t',index=False)
    build_pdf(markers,counts,q)
    (OUT/'source_checksums.json').write_text(json.dumps(SOURCES,indent=2)+'\n')
    print(PDF)
    print('Five-page review generated; all 3,434 MG cell IDs checked for uniqueness and sample/donor/eye consistency.')


if __name__=='__main__':main()
