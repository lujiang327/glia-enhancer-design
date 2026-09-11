#!/usr/bin/env python3
"""Build the feasibility report from measured outputs; fail on missing inputs."""
import csv
import json
from pathlib import Path
import statistics
import base64
import html
import re


def table(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))


def run():
    q=Path('reports/qc')
    audit=json.loads((q/'fragment_audit.json').read_text());pool=json.loads((q/'pool_qc.json').read_text())
    donors=table(q/'donor_summary.tsv');frip={r['scope']:r for r in table(q/'frip.tsv')}
    peak=json.loads((q/'peak_qc.json').read_text());corr=table(q/'donor_correlations.tsv')
    offdiag=[float(r['pearson_log1p_CPM_50kb']) for r in corr if r['donor_a']!=r['donor_b']]
    n=sum(int(r['unambiguous_MG_cells']) for r in donors);depth=sum(int(r['retained_fragment_records']) for r in donors)
    unique={r['bare_barcode'] for r in table(q/'barcode_audit.tsv') if r['n_matching_RNA_samples']=='1'}
    retained={r['bare_barcode']:int(r['retained_records']) for r in table(q/'filtered_barcode_counts.tsv')}
    med=statistics.median(retained.get(b,0) for b in unique)
    cohort_frip=sum(int(frip[r['donor']]['fragments_overlapping_peak']) for r in donors)/depth
    t='| Donor | MG cells* | Retained fragments* | Median/cell* | TSS proxy* | FRiP* |\n|---|---:|---:|---:|---:|---:|\n'
    for r in donors:
        d=r['donor'];t+=f"| {d} | {int(r['unambiguous_MG_cells']):,} | {int(r['retained_fragment_records']):,} | {float(r['median_retained_fragments_per_cell']):,.0f} | {float(r['pooled_TSS_proxy']):.2f} | {float(frip[d]['FRiP']):.1%} |\n"
    report=f'''# Müller glia ChromBPNet feasibility — GSE196235

Assessed 2026-09-09. **GO for a training pilot using the donor-resolved subset. NO-GO for differential accessibility, design, or library preparation until the model is trained and validated.** Current model status: NOT TRAINED. This recommendation is a data-feasibility judgment, not proof of a reliable model.

## Measured depth and quality

The working cohort contains **{n:,} published-label Müller glia cells**, **{depth:,} retained fragment records**, a median of **{med:,.1f} fragments/cell**, and cohort FRiP **{cohort_frip:.1%}**. Four donor pseudobulks are saved. Both eyes are combined within biological donor.

{t}
*All donor values describe the **unambiguously mapped subset**, not complete donor totals. Complete Müller glia cell counts and depths per donor remain unresolved for excluded barcodes; they are not zero. “Retained” means chr1–22,X, valid coordinates, no overlap with hg38 blacklist v2. Original alignment quality/deduplication is inherited from the processed release, not revalidated from BAMs.*

The complete published MG pool contains {audit['fragment_records']:,} records; {pool['counts']['retained']:,} survive genomic filters. Read-support sum is {audit['read_support_sum']:,}; it is **not** independent fragment depth. There are {audit['distinct_bare_barcodes']:,} distinct bare barcodes, which cannot be equated to biological cells. Full-pool TSS proxy is {pool['pooled_tss_enrichment_proxy']:.2f}; full-pool FRiP is {float(frip['full_pool']['FRiP']):.1%}.

MACS3 called **{peak['diagnostic_peaks_after_blacklist']:,} diagnostic accessible peaks** after blacklist filtering ({peak['diagnostic_peaks_before_blacklist']:,} before). These are full-pool BEDPE peaks at q=0.01, not DA regions or the final relaxed ChromBPNet training peak set. FRiP counts each retained fragment once if it overlaps any diagnostic peak. The TSS proxy averages +/-50 bp relative to outer 100-bp flanks in +/-2 kb around {pool['tss_count']:,} unique strand-specific RefGene TSSs; this differs from ArchR per-cell TSS scores.

![Donor counts, depth, TSS profiles and consistency](figures/donor_qc.png)

## Annotation, donor consistency and quality caveats

Published labels are retained through membership in the authors' MG fragment file. Matching paired RNA supports them: GLUL is detected in 97.1% and RLBP1 in 86.7% of the working cohort. GLUL/RLBP1 mean normalized expression is higher than the other-barcode background. SLC1A3, SOX9, AQP4 and retinal off-target markers are included in the saved marker table. No cells were manually relabeled.

Bare barcodes lack sample prefixes: **231 barcodes and 6,288,321 records (6.76% of the full pool)** match multiple RNA samples. Excluding them preserves {n:,} identifiable cells and avoids guessing donors. The archived authors' RNA code supplies prefix-to-sample mapping; GEO/HCA supplies donor metadata. The four demographics are LVG1: 55/M; LGS1: 85/F; LGS2: 87/M; LGS3: 74/F with dementia recorded, although retinal tissue is marked normal. [GSE196235](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196235), [HCA project](https://explore.data.humancellatlas.org/projects/4f4f0193-ede8-4a82-8cb0-7a0a22f06e63), [archived annotation code](https://doi.org/10.5281/zenodo.6795162).

Donor ATAC correlations are {min(offdiag):.3f}–{max(offdiag):.3f} (Pearson, log1p CPM in blacklist-free 50-kb midpoint bins with nonzero pooled counts). TSS profiles agree and fragment lengths show nucleosomal structure. These coarse correlations do **not** establish base-resolution generalization or absence of batch effects. Donor, age and technical effects are confounded.

Neuronal RNA is a material caveat: PDE6A is detected in 55.5% overall, ranging from 36.7% in LVG1 to 75.6% in LGS3. Ambient RNA and residual mixed nuclei are possibilities, not established explanations. Require held-out-donor testing, rod/glial locus and motif inspection, and sensitivity to excluding LGS3 before calling the model reliable. Do not infer a new label from these markers alone.

![RNA marker support and fragment lengths](figures/markers_and_lengths.png)

## Available files and fallback

GSE196235 provides hg38 paired multiome data, 13 cell-type fragment pools and a 36,601-feature × 51,645-cell RNA matrix. The MG file is physically attached to GSM5866073, despite series prose pointing to rep1; it pools replicates. Downloaded MG fragments and all three RNA matrix components match HCA SHA256 checksums. No per-cell ATAC QC table or peak-by-cell matrix was present in the inspected GEO supplementary inventory. Therefore original per-cell TSS/FRiP distributions and full donor totals are unavailable from these inspected assets. Published filtering included TSS >6 and >2,500 ATAC fragments; that is a reported selection rule, not a newly measured per-cell result. [ATAC sample methods](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5866079).

The 2022 release includes **legacy BPNet** code/models and MG predictions, a useful precedent but not current ChromBPNet validation. [Study](https://doi.org/10.1016/j.xgen.2022.100164), [model code](https://doi.org/10.5281/zenodo.6796067).

Additional depth is not the current limiting factor. GSE265801 is a heterogeneous SuperSeries; do not pool it wholesale. Its GSE265774 archive listing contains 63 matrix files and no fragment-named files. Adult subseries GSE281526 has 31 fragment files totaling 57.07 GB; GSM8622730 explicitly supplies GRCh38/GENCODE32 multiome fragments, index and donor BCM_22_0047 metadata. Older GSM7064216 is hg19. The queried HRCA CELLxGENE collection contains RNA datasets, which cannot replace ATAC fragments. Atlas cell labels still need a validated barcode join before expansion. Avoid counting reanalyzed GSE196235 nuclei twice; exclude fetal cohorts from an adult pilot. [SuperSeries](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE265801), [adult fallback sample](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8622730), [older hg19 sample](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7064216), [atlas resources](https://rchenlab.github.io/resources/human-atlas.html).

## Next gate and execution blocker

Proceed to a genome-wide ChromBPNet pilot with the donor-resolved cohort, matched genomic nonpeaks and the saved chromosome folds. First verify the reference/coordinate convention, generate insertion tracks, finalize relaxed training peaks/backgrounds, and fit a multiome-compatible bias model. Then assess multiple seeds, held-out counts/profiles, donor generalization and attribution/motif recovery. The current data do not show an obvious depth-based reason to reject a pilot. [Official preprocessing](https://github.com/kundajelab/chrombpnet/wiki/Preprocessing), [bias-model guidance](https://github.com/kundajelab/chrombpnet/wiki/Bias-model-training).

This workspace is macOS ARM64 with no detected NVIDIA tooling or ChromBPNet installation. No Linux/GPU execution target has been supplied. Training, bias fitting, validation losses, profile/count prediction metrics, attributions and checkpoints are therefore **not available**. No differential accessibility, candidate ranking, sequence optimization or brain-astrocyte analysis has run. See [workflow and validation criteria](../docs_workflow.md), [training parameters](../config/training_plan.json), and [validation gate](../config/model_validation.json).
'''
    Path('reports/phase1_qc_report.md').write_text(report)
    # Self-contained HTML for this report's deliberately small Markdown subset.
    def inline(text):
        text=html.escape(text)
        text=re.sub(r'\*\*(.+?)\*\*',r'<strong>\1</strong>',text)
        text=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',text)
        return text
    blocks=[]
    for block in report.strip().split('\n\n'):
        if block.startswith('# '):blocks.append('<h1>'+inline(block[2:])+'</h1>')
        elif block.startswith('## '):blocks.append('<h2>'+inline(block[3:])+'</h2>')
        elif block.startswith('|'):
            lines=block.splitlines();head=[v.strip() for v in lines[0].strip('|').split('|')]
            body=''.join('<tr>'+''.join('<td>'+inline(v.strip())+'</td>' for v in l.strip('|').split('|'))+'</tr>' for l in lines[2:])
            blocks.append('<table><thead><tr>'+''.join('<th>'+inline(v)+'</th>' for v in head)+'</tr></thead><tbody>'+body+'</tbody></table>')
        elif block.startswith('!['):
            m=re.fullmatch(r'!\[([^\]]*)\]\(([^)]+)\)',block)
            data=base64.b64encode((Path('reports')/m[2]).read_bytes()).decode()
            blocks.append('<figure><img alt="'+html.escape(m[1])+'" src="data:image/png;base64,'+data+'"></figure>')
        else:blocks.append('<p>'+inline(block)+'</p>')
    page='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Müller glia QC and model feasibility</title><style>body{max-width:1050px;margin:40px auto;padding:0 24px;font:16px/1.6 system-ui;color:#19313b}h1{font-size:30px;line-height:1.2}h2{margin-top:32px;color:#246b68}table{width:100%;border-collapse:collapse;font-size:14px}th,td{padding:9px;border-bottom:1px solid #dce3e6;text-align:right}th:first-child,td:first-child{text-align:left}figure{margin:24px 0}img{width:100%;height:auto}a{color:#176c95}@media print{body{margin:0;font-size:11px}figure,table{break-inside:avoid}}</style><body>'+''.join(blocks)+'</body></html>'
    Path('reports/phase1_qc_report.html').write_text(page)
    status=dict(data_recommendation='GO_PILOT',model_reliability='NOT_ASSESSED_NOT_TRAINED',downstream_recommendation='NO_GO',working_cells=n,working_fragment_records=depth,working_median_fragments=med,working_FRiP=cohort_frip,execution_blocker='Linux/GPU training environment not configured',caveats=['ambiguous barcodes excluded','neuronal RNA and donor-specific variation','LGS3 dementia metadata','QC proxies differ from original per-cell metrics'])
    Path('config/qc_recommendation.json').write_text(json.dumps(status,indent=2)+'\n')
    print(json.dumps(status,indent=2))


if __name__=='__main__':run()
