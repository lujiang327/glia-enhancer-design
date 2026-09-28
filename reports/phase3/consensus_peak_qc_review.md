# Phase 3 consensus-peak QC review

## Decision

**Accepted for donor-level counting. Differential accessibility remains
blocked pending common-universe count-matrix QC.**

The corrected build contains 177,777 nonoverlapping 500-bp peaks selected from
582,543 pooled candidates that were supported by at least two donor peak calls
in the same published cell type. The build completed as Great Lakes job
61889009 from commit `0a5eb02`; all four output checksums passed.

Independent checks found zero malformed widths, overlapping or duplicate
intervals, rank quantiles outside `[0,1]`, and hg38-blacklist overlaps. Donor
support is:

| Supporting donors | Peaks | Fraction |
|---:|---:|---:|
| 2 | 45,316 | 25.49% |
| 3 | 36,243 | 20.39% |
| 4 | 96,218 | 54.12% |

## Cross-cell-type center selection

Raw MACS3 q-values scale with cell-type depth. The final build therefore ranks
q-value and signal within each cell type and donor-support stratum before
resolving overlaps between cell types. Raw q-value remains only a deterministic
tie-breaker. Müller glia supplies 42,086 final peak centers (23.67%) and rod
supplies 65,698 (36.95%). In the archived raw-q build, these values were 37,971
(21.36%) and 75,541 (42.49%), respectively. The correction reduces the depth
advantage of the rod pool while retaining donor support as the first priority.

The corrected and archived universes cover the same loci at this resolution:
all 177,777 corrected peaks overlap an archived peak, and all 177,776 archived
peaks overlap a corrected peak. The correction changes which summit represents
an overlapping locus rather than introducing a new set of biological regions.

## Remaining gate

The accepted universe is suitable for constructing the 52-column donor-by-cell-
type count matrix. The next review must examine FRiP, detected peaks, donor
correlations within each published cell type, PCA cell-type versus donor effects,
the low-cell astrocyte/microglia/RGC strata, and rod depth imbalance. No
differential model may be fit before that review is complete.
