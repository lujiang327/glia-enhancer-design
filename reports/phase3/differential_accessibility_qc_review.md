# Phase 3 donor-blocked differential-accessibility review

## Decision

**GO for ranking approximately 1,000 Muller-glia-accessible regulatory
regions, with donor and glial-specificity safeguards.** This decision opens
candidate ranking and annotation. It does not select candidates, establish
enhancer activity, or authorize sequence optimization before the ranked set is
reviewed.

Great Lakes edgeR job `62403767` completed with exit code 0 from commit
`0813049`. All 13 compact output checksums pass. The complete contrast and
leave-one-donor-out result tables remain under
`PHASE3_SCRATCH_ROOT/differential_accessibility/`; the compact evidence is in
`reports/phase3/differential_accessibility_qc/`.

## Model and normalization QC

The analysis retained 147,359 of 177,777 consensus peaks (82.89%) after
design-aware `filterByExpr`. The donor-blocked design has 52 samples, 16
coefficients, rank 16, and condition number 17.56. The design is estimable and
does not show a rank or severe-conditioning failure.

TMM normalization factors range from 0.661 to 1.419 and have product 1. The
largest factor belongs to the low-depth RGC LGS1 sample and the smallest to the
very deep rod LGS2 sample. The factors respond to known depth and composition
differences without an isolated unexplained extreme. All 52 samples remain in
the analysis.

The common dispersion is 0.0426 (BCV 0.206). The tagwise BCV and quasi-
likelihood dispersion plots show smooth abundance-dependent trends with robust
shrinkage and a limited high-dispersion tail. MDS separates the published
retinal cell types while keeping biological replicates broadly coherent. No
global donor axis or failed sample is evident.

P-value histograms contain strong near-zero enrichment for the neuronal
contrasts, as expected for distinct retinal cell types. The closer astrocyte
and microglia comparisons retain broad null-like tails plus a near-zero signal
component. There is no general excess at intermediate small p-values suggesting
an obvious calibration failure. The primary contrast also has a conservative
mass near one; this does not create false-positive enrichment, but candidates
must be ranked by effect size and donor evidence as well as FDR.

## Differential signal

The equal-cell-type-weight Muller-versus-rest contrast reports 79,471 peaks at
FDR below 0.05: 51,821 are Muller-positive and 27,650 Muller-negative. Of the
positive set, 34,030 have log2 fold change at least 1. Separately, 47,496
significant positive peaks have a positive within-donor target-minus-rest
effect in all four donors.

The closest glial comparisons contain fewer Muller-positive peaks than the
neuronal comparisons, as biologically expected:

| Contrast | Positive, FDR < 0.05 | Positive and log2FC >= 1 | Positive in all donors |
|---|---:|---:|---:|
| Muller vs astrocyte | 11,612 | 9,384 | 10,926 |
| Muller vs microglia | 16,718 | 14,770 | 16,335 |

All twelve prespecified pairwise comparisons completed. Rod depth does not
determine the primary result because each off-target cell type has equal weight
in the primary contrast.

## Donor sensitivity

For Muller versus astrocyte, leave-one-donor-out logFC correlations with the
four-donor fit range from 0.810 to 0.954; for Muller versus microglia they range
from 0.832 to 0.960. Every full-model Muller-positive peak retains a positive
effect direction in every leave-one-donor-out refit.

FDR retention is less stable than effect direction. It ranges from 32.2% to
99.9% for astrocytes and 59.2% to 97.8% for microglia. The weakest cases occur
after removing LVG1 from the astrocyte comparison and LGS3 from the microglia
comparison. With only three donors remaining, this loss can reflect reduced
power as well as donor influence. Candidate ranking must therefore retain the
minimum leave-one-donor-out logFC and maximum leave-one-donor-out FDR for both
glial comparisons. Stable positive direction is required; universal
leave-one-donor-out FDR significance is a ranking preference rather than a hard
filter.

## Candidate-ranking contract

Candidate ranking may proceed using the full 147,359-peak tested universe. The
primary eligibility screen is:

1. Muller-versus-rest FDR below 0.05 and log2FC at least 1;
2. positive within-donor Muller-versus-rest effect in all four donors;
3. positive full-model log2FC against astrocyte and microglia; and
4. positive leave-one-donor-out log2FC in all eight glial sensitivity fits.

Eligible regions should then be ranked using primary effect size and FDR,
minimum effect across all 12 pairwise comparisons, target accessibility,
maximum off-target accessibility, glial leave-one-donor-out stability, and
donor reproducibility. The final approximately 1,000-region set must preserve
genomic-class annotations and should not be described as confirmed enhancers.

Regions that are strongly Muller-versus-rest but accessible in astrocytes or
microglia may remain in the full scored table, but they should not outrank
regions with comparable target activity and stronger glial specificity. Keep
the full scored universe and all component metrics so the ranking can be
audited and thresholds can be changed without rerunning edgeR.

## Remaining caveats

- Four donors limit power and prevent separating all age, donor, and technical
  effects.
- LGS3 has dementia metadata despite the retina being marked normal.
- Low-cell astrocyte, microglia, and RGC strata increase uncertainty in the
  corresponding pairwise comparisons.
- Leave-one-donor-out analysis refits the differential model; it does not
  retrain ChromBPNet without that donor.
- The Muller ChromBPNet model retains a nonpeak count-calibration caveat.
- No validated off-target retinal ChromBPNet models are currently available;
  observed off-target accessibility must initially provide the specificity
  penalty.

