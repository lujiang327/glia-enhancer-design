# GSE196235 processing, clustering and annotation review

Prepared October 20, 2026 for discussion with Thanh. The five-page handout is
`output/pdf/MG_multiome_processing_annotation_review.pdf` from the repository
root. This review documents existing analyses, inspects archived code, checks
retained MG cell-ID mappings, and replots saved QC results. It does not rerun
alignment, clustering, doublet detection, or differential analysis.

## Main conclusion

We used the original study's cell-type annotations. We did not independently
recreate the final single-cell clustering or relabel cells. Published markers,
our matched-RNA MG check, and donor-level ATAC diagnostics provide support, but
they cannot certify that all upstream processing or individual labels are
issue-free. Unresolved neuronal RNA in MG-labelled nuclei, limited off-target
glial depth, and incomplete independent annotation auditing remain material.

## Original study: reported methods

Wang et al., *Cell Genomics* (2022),
[doi:10.1016/j.xgen.2022.100164](https://doi.org/10.1016/j.xgen.2022.100164),
profiled eight postmortem retinas from four individuals with joint 10x multiome.
The final published RNA dataset has 51,645 nuclei in 22 clusters, grouped into
13 cell types. The figure in this packet is their published Figure 1, not a new
embedding from our project.

| Step | Published method / parameter | Our replication status |
|---|---|---|
| Alignment and count generation | Cell Ranger ARC 2.0.0, hg38 | Inherited processed release; no FASTQ/BAM reanalysis |
| Paired-nucleus selection | Barcodes present in RNA and ATAC | Inherited original paired-cell filtering |
| Joint QC | RNA UMIs >200 and <50,000; mitochondrial ratio <1%; ribosomal ratio <5%; ATAC TSS >6; ATAC fragments >2,500 | Published thresholds, not newly measured per-cell distributions |
| Automated doublets | ArchR addDoubletScores / filterDoublets | Inherited; not rerun |
| RNA preprocessing | Seurat 3.1.5; NormalizeData; 2,000 variable genes; ScaleData; PCA | Archived script inspected; not rerun |
| Sample correction | Harmony using eight sample IDs | Archived initial clustering explicitly uses Harmony |
| RNA clustering | Graph neighbors on first 20 dimensions; resolution 0.5 | Original final object / exact execution not reproduced |
| Cluster review | Marker-based annotation; remove mixed-marker and marker-poor clusters, then recluster | Authors' decisions inherited; no manual relabeling by us |
| Final labels | 22 clusters grouped into 13 types | Retained through membership in published fragment files |
| RNA-to-ATAC annotation | Paired-nucleus barcodes; ArchR receives RNA identities | Original transfer not independently replayed |

The inspected workflow does not establish that our project performed joint WNN
clustering or atlas-based label transfer. In the source study, RNA identities
were carried to paired ATAC nuclei. Its later ATAC LSI/UMAP is a separate
accessibility embedding.

### Published marker examples

These are examples in the authors' RNA code / figure, not exclusive marker rules
or a newly fitted classifier.

| Published type | Example genes |
|---|---|
| Muller glia | GLUL, RLBP1 |
| Astrocyte | PAX2, GFAP |
| Microglia | C1QA, HLA-DRA |
| Rod | PDE6A, NR2E3 |
| Cone | ARR3, OPN1LW |
| ON cone bipolar | ISL1, GRM6 |
| OFF cone bipolar | GRIK1 |
| Rod bipolar | PRKCA |
| GABA amacrine | GAD1, GAD2 |
| Glycinergic amacrine | SLC6A9 |
| AII amacrine | CALB2, GJD2 |
| Horizontal | ONECUT1, ONECUT2 |
| Retinal ganglion | NEFL, NEFM |

## Our processing of the released files

1. **Input provenance.** Use GSE196235 processed hg38 data: the 36,601-feature by
   51,645-cell RNA matrix and 13 cell-type ATAC fragment pools. The pooled RNA
   components are attached to GSM5866081; cell-type fragment pools are attached
   to GSM5866073. Those attachment names must not be mistaken for a one-eye or
   one-donor cohort. The files contain data across the eight samples. We did
   not add GSE265801 or brain-astrocyte data to this cohort.
2. **Sample and donor recovery.** RNA IDs contain sample prefixes 1-8, whose
   mapping is explicit in archived author code. Combine with public GEO/HCA
   donor metadata. Prefixes 1/2 are LVG1, 3/4 LGS1, 5/6 LGS2, 7/8 LGS3;
   odd/even entries correspond to the two eyes. See `config/samples.tsv`.
3. **Resolve only unique matches.** Published ATAC pools have bare 10x barcodes
   without sample prefixes. Compare against all prefixed RNA IDs. Assign a
   donor only when a bare barcode matches exactly one RNA sample. Do not guess
   a donor when a barcode is reused. For MG, exclude 231 ambiguous bare barcodes
   representing 6,288,321 fragment records; retain 3,434 uniquely resolved MG
   nuclei. The complete original donor totals for excluded barcodes remain
   unresolved. Distinct bare barcodes are not complete biological cell counts.
4. **Retain published type membership.** Membership in a cell-type fragment file
   supplies the published label. Our scripts did not create RNA clusters,
   transfer inferred labels, or change cell identities. All 13 types have
   retained cells from four donors; the existing audit totals 48,422
   donor-assigned cell records across the pools. A fresh global cross-pool
   conflict check is still needed to audit those memberships jointly.
5. **RNA MG support.** Match the unambiguous MG IDs to RNA matrix columns and
   compute detection and mean log1p(CP10k) for the 12 saved genes. This is a
   targeted MG marker check against a mixed background, not an independent
   all-13-type marker audit or clustering validation.
6. **ATAC filtering and depth.** Keep chr1-22/X, valid coordinates, and fragments
   not overlapping hg38 blacklist v2. Count each retained fragment row once;
   do not expand fifth-column read support into independent observations.
   Processed alignment and deduplication are inherited. Both eyes are combined
   within each biological donor.
7. **Preserve biological replicates.** Create 52 donor-by-type pseudobulks. Total
   retained fragments are 775,669,056, including 85,979,585 MG fragments.
   Review fragment lengths, pooled TSS proxy, FRiP, depth, correlations and PCA.
   Our TSS proxy is not the original per-cell ArchR statistic. Diagnostic MG
   FRiP uses a different peak universe from later consensus-matrix FRiP.
8. **Downstream order.** Train and validate MG ChromBPNet first using genome-wide
   regions. Then conduct donor-blocked DA over a common 177,777-region consensus
   universe constructed across all 13 types, with 147,359 regions retained for
   testing. The model was not trained only on the eventual 1,000 candidates.
   Existing DA uses TMM and donor blocking, with an equal-cell-type-weight
   primary contrast and glial donor-omission sensitivity checks.

## Checks completed in this review

- The 3,434 retained MG `cell_id` values are unique, as are their bare barcodes.
- Each retained MG ID exactly matches its sample prefix plus bare barcode.
- Every MG record agrees with the archived sample-prefix-to-donor and eye map.
- The 52-row count table sums to 48,422 donor-assigned cell records; MG sums to
  3,434. These agree with existing reports. This arithmetic does not replace
  a global cross-pool identity check.
- All five final PDF pages were rendered and visually inspected. Saved figure
  and PDF inputs have SHA256 provenance.

These checks validate internal identifier consistency, not the original cell
biology or every step of the authors' label transfer. A source prefix error or
an inherited label error could remain internally consistent.

## Known caveats and exact code references

### Rod RNA in MG-labelled cells

Our existing marker summary detects GLUL in 97.1% and RLBP1 in 86.7% of retained
MG nuclei, supporting MG identity. It also detects PDE6A in 55.5% and RHO in
63.7%. LGS3 shows the strongest PDE6A detection (75.6%). Ambient RNA and residual
mixed nuclei are possibilities, not demonstrated explanations. No ambient-RNA
correction or new doublet inference was performed. Detection of a rod transcript
alone is not sufficient grounds to relabel or discard a nucleus. The next
assessment should consider coexpression, relative counts and donor dependence.

### Sparse glial off-target strata and donor confounding

Astrocyte LGS1/LGS2 contain 12/46 cells, microglia LVG1 contains 27, and RGC LGS1
contains 42. Absent browser signal is weak evidence of biological absence in
these strata. Broad ATAC sample structure is coherent, but some glial samples
correlate most strongly with same-donor MG. Donor blocking and sensitivity
analyses help address this; they do not prove perfect annotation. LGS3 also has
dementia recorded despite normal-retina metadata. Age, donor, eye and technical
effects cannot all be independently separated with this cohort.

### Archived-code reproducibility questions, not demonstrated annotation errors

The archived code is under Zenodo
[doi:10.5281/zenodo.6795162](https://doi.org/10.5281/zenodo.6795162), extracted
locally to `data/metadata/author_code/seankwang-retina-multiome-651e181/`.

- `scRNA:143-162`: initial Harmony correction and neighbor/UMAP construction
  explicitly specify the Harmony reduction.
- `scRNA:184-192`: manual cluster exclusion and a later reclustering block omit
  an explicit reduction argument. The precise final execution requires the
  final object/logs; the initial Harmony call alone cannot establish which
  reduction was used for every subsequent operation.
- `scRNA:200-203`: the ordered 22-entry identity vector is applied to the
  existing cluster levels. We have not reconstructed those cluster levels.
- `scRNA:235-243`: prefix 1-8 to sample-name mapping is explicit.
- `scATAC:75-78`: RNA clusters and identities are assigned as positional vectors
  to the ArchR project. This snippet has no visible named-cell reorder
  assertion. One-to-one barcode identity and order should be checked against
  the final objects before independently replaying it.

We did not execute those snippets as a fresh annotation pipeline; we used the
published labelled fragment pools. These observations limit reproducibility
claims but do not demonstrate that the authors' released annotations are wrong.

## Most useful next annotation audit

1. Obtain the sample-prefixed per-cell identities, final RNA object/embedding,
   and per-cell ATAC QC/doublet information where available.
2. Check global uniqueness across all 13 cell pools, explicit RNA/ATAC cell-ID
   joins, sample prefix, donor and eye, with no positional label assignment.
3. Reproduce a canonical-marker dot plot across all 13 types, plus MG versus
   astrocyte/microglia views separately by donor. Distinguish a reproduced
   embedding from a new clustering run.
4. Inspect rod-high MG nuclei for coexpression and RNA count burden. If the
   evidence supports a contamination subset, quantify sensitivity of existing
   pseudobulks and selected parents before changing any labels or exclusions.
5. Retrieve/reproduce original QC and doublet diagnostics if a stronger claim
   about raw processing is needed. Successful model fitting or DA is not a
   substitute for these checks.

## Files and reproduction

- `01_processing_workflow.png`: published versus our workflow.
- `02_published_Wang2022_Figure1.jpg`: unchanged original UMAP/marker figure.
- `03_MG_RNA_marker_review.png`, `RNA_marker_evidence.tsv`: our existing MG marker
  QC replotted without per-gene color scaling.
- `04_donor_celltype_ATAC_review.png`, `donor_celltype_counts.tsv`: donor coverage
  and ATAC pseudobulk PCA, not a single-cell embedding.
- `RNA_and_MG_sample_counts.tsv`: RNA-column and retained-MG counts by sample.
- `validation.json`, `source_checksums.json`, `output_checksums.sha256`: checks
  and reproducibility records.

```bash
MPLCONFIGDIR=/tmp/glia-matplotlib .venv/bin/python scripts/build_annotation_process_review.py
```

The source Figure 1 was downloaded from the public PMC CDN and is reproduced
unchanged with attribution under CC BY 4.0. Source:
https://pmc.ncbi.nlm.nih.gov/articles/PMC9584034/#fig1.
It is published contextual evidence, not a figure generated from our retained
cohort. The original PDF-quality source may provide sharper lettering than the
748-pixel public JPG used here.
