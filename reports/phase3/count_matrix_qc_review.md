# Phase 3 common-universe count-matrix QC review

## Decision

**GO for donor-blocked differential accessibility with caveats.**

The 52 donor-by-published-cell-type samples were counted on the accepted
177,777-peak retinal universe. Great Lakes array job `62258579` and dependent
summary job `62258580` completed with exit code 0. All 52 sample checksum sets,
the full matrix checksum, and all compact report checksums passed. The full raw
count matrix remains in scratch and the compact QC outputs are versioned under
`reports/phase3/count_matrix_qc/`.

No differential test was run during this QC stage.

## Depth and peak occupancy

- Total retained fragment rows: 775,669,056.
- Fragments overlapping at least one consensus peak: 402,711,528 (weighted
  FRiP 0.5192).
- Per-sample FRiP range: 0.2988–0.6328.
- Per-sample detected-peak range: 69,281–174,869.
- Peaks passing the exploratory filter of CPM >= 1 in at least four samples:
  177,342 of 177,777 (99.76%). This filter was used only for correlation and
  PCA QC; differential filtering must be recalculated with the edgeR design.
- Fragments intersecting more than one nonoverlapping peak: 1.379% of retained
  rows. These produced 10,725,535 additional peak assignments and are retained
  under the documented interval-overlap counting policy.

Müller-glia donor pseudobulks contain 16.3–29.4 million retained fragments,
have FRiP 0.5466–0.5609, and detect 160,668–172,169 peaks. Their within-cell-type
donor correlations range from 0.8652 to 0.9327 (median 0.9009).

## Biological and donor structure

PCA on the 20,000 most variable exploratory log2-CPM peaks separates the
published retinal cell types. PC1 and PC2 explain 45.17% and 20.95% of variance.
Cell-type eta-squared is 0.9705 and 0.9727 for these components, whereas donor
eta-squared is 0.0049 and 0.0006. Cell type remains the dominant association
through PC8; there is no evidence that a global donor batch dominates the
matrix.

Within-cell-type median donor correlation ranges from 0.7058 for microglia and
0.7114 for AII amacrine to 0.9594 for rods. Six of 52 samples have a highest
correlation with a biologically related different group rather than their own
published group: AII amacrine LGS3 with GABA amacrine LGS3; astrocytes LGS2 and
LGS3 with Müller glia from the same donor; and microglia LGS1–LGS3 with Müller
glia from the same donor. The PCA still separates these cell types. These
relationships support donor blocking and careful pairwise review rather than
sample removal or manual relabeling.

Within the three glial groups, same-donor correlations are modestly higher than
cross-donor correlations: median differences are 0.0326 for Müller glia versus
astrocyte, 0.0645 for Müller glia versus microglia, and 0.0363 for astrocyte
versus microglia. Donor blocking is therefore essential, and the Müller-glia
versus astrocyte and Müller-glia versus microglia contrasts require leave-one-
donor-out sensitivity checks before they can support specificity claims.

## Samples retained with caveats

No sample is excluded. Astrocyte LGS1 is the smallest pseudobulk (434,120
fragments) but has FRiP 0.4921 and 69,281 detected peaks. RGC LGS1 has the lowest
FRiP (0.2988) but 1,760,704 fragments and 108,478 detected peaks; the RGC median
donor correlation is 0.7806. The lower-depth astrocyte, microglia, RGC, AII
amacrine, and glycinergic amacrine strata must remain visible in diagnostics.

Rod depth reaches 144.3 million fragments in LGS1 and is far above the other
cell types. The primary Müller-glia-versus-rest contrast must give cell types
equal biological weight rather than letting rod abundance determine the
contrast.

## Requirements for differential accessibility

The next stage may fit donor-blocked edgeR quasi-likelihood models using raw
counts. It must:

1. retain all four biological donors and all 52 pseudobulks;
2. use a design with donor and published cell type;
3. calculate sample-specific edgeR normalization factors and inspect their
   relationship with library depth and cell type;
4. use design-aware abundance filtering rather than the exploratory QC filter;
5. fit the primary Müller-glia versus equal-weight average of the 12 off-target
   cell types and the 12 prespecified pairwise contrasts;
6. use robust dispersion estimation and report model/design rank, BCV,
   mean-variance behavior, p-value calibration, and donor consistency;
7. preserve low-depth samples and report sensitivity results if any sample has
   disproportionate influence; and
8. run leave-one-donor-out sensitivity checks for the Müller-glia versus
   astrocyte and Müller-glia versus microglia contrasts; and
9. stop before candidate selection if normalization, dispersion, or donor
   diagnostics are unacceptable.
