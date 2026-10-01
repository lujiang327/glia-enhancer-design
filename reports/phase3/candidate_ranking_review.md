# Muller-glia candidate-ranking review

## Decision

**GO for ChromBPNet parent-sequence scoring of the provisional top 1,000
regions, with caveats.** Candidate ranking job `62978649` completed with exit
code 0 from commit `c20c16e` in 51 seconds using 365 MB peak memory. All four
Great Lakes output checksums passed. Variant generation remains blocked until
the parent sequences and model predictions pass their own QC.

The full scored universe remains in scratch as
`candidate_ranking/all_scored_regions.tsv.gz` with SHA256
`bb46032fa97bc327ee446185ea03699272af1c7c2a47415e7727524f59f9c844`.
The compact summary, top table, BED file, and checksums are versioned under
`reports/phase3/candidate_ranking/`.

## Eligibility and ranking

Of 147,359 tested regions, 21,269 passed the prespecified eligibility gate:

- Muller-versus-rest FDR below 0.05 and log2FC at least 1;
- positive Muller-versus-rest effect in all four donors;
- positive full-model effects against astrocytes and microglia; and
- positive effects in all eight glial leave-one-donor-out refits.

The top 1,000 are unique, nonoverlapping 500-bp regions. Ranks are complete and
scores decrease monotonically. Every selected region is positive in all 12
full-model pairwise contrasts; this emerged from the ranking rather than being
an extra hard filter. Eight hundred seventy-two also remain significant at FDR
below 0.05 in every glial leave-one-donor-out refit. All 1,000 retain positive
leave-one-donor-out effect direction by construction.

Selected-region metric ranges are:

| Metric | Minimum | Median | Maximum |
|---|---:|---:|---:|
| Primary Muller-vs-rest log2FC | 3.399 | 5.205 | 6.883 |
| Minimum pairwise log2FC | 1.377 | 1.955 | 3.741 |
| Minimum four-donor primary effect | 2.767 | 3.791 | 5.102 |
| Minimum glial leave-one-donor-out log2FC | 1.033 | 1.669 | 3.786 |
| Observed accessibility specificity, log2 ratio | 1.294 | 1.938 | 3.625 |
| Muller mean CPM | 6.443 | 18.424 | 68.655 |
| Maximum off-target mean CPM | 0.733 | 4.626 | 20.364 |

The largest observed off-target is microglia for 619 regions and astrocytes for
301, consistent with the need to prioritize specificity among related glial
cell types. The remaining 80 regions have a neuronal cell type as their largest
off-target.

## Genomic composition

The top set retains the intended mixture of open regulatory-region classes:

| Genomic class | Regions |
|---|---:|
| Distal intergenic | 501 |
| Intronic | 420 |
| Exonic | 43 |
| Promoter | 36 |

All canonical chromosomes are represented. No region is labeled a confirmed
enhancer. RefGene class is a coordinate annotation with promoter priority; it
does not establish regulatory activity or target-gene assignment. Nearest-gene
labels are descriptive and must not be treated as validated enhancer-gene
links.

## Phase 4 requirements

The next step may extract the natural hg38 sequence for every selected region
and score it with the frozen validated Muller no-bias ChromBPNet checkpoint.
The parent table must retain candidate rank, peak ID, coordinates, sequence,
differential statistics, genomic class, observed specificity, and checksum.

Before designing variants, review sequence extraction, reference allele and
length integrity, invalid-base frequency, prediction distributions, edge
effects, and agreement between predicted accessibility and observed Muller
accessibility. Sequence attribution and perturbation should use the same frozen
checkpoint and preserve the unmodified parent for every future variant.

## Caveats

- The composite score and its weights are prespecified prioritization choices,
  not a statistically calibrated probability of enhancer activity.
- The top set emphasizes strong cross-cell-type specificity; canonical Muller
  marker loci need not rank highly if they are accessible in related glia.
- Off-target specificity currently uses observed pseudobulk accessibility
  because validated off-target ChromBPNet models are unavailable.
- Four donors limit sensitivity analysis, especially for astrocytes and
  microglia.
- ChromBPNet retains the previously documented nonpeak count-calibration
  caveat.

