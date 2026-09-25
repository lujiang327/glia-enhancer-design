# Phase 3 donor pseudobulk QC review

**Decision: GO for reproducible peak calling, with caveats. Differential
accessibility remains blocked until the common peak universe and its count-matrix
QC pass review.**

All 52 expected donor-by-cell-type pseudobulks completed successfully. Together
they contain 775,669,056 retained fragment rows and 1,163,781,495 fifth-column
read-support counts after restriction to chr1-22/X, coordinate validation, and
full-fragment exclusion against the hg38 blacklist. The smallest pseudobulk is
astrocyte LGS1 with 434,120 retained fragment rows; no sample fell below the
prespecified 250,000-fragment screen.

Blacklist loss is low and consistent: 5,667,245 assigned canonical fragment rows
were excluded, or 0.725% overall. Cell-type-specific loss ranges from 0.639% for
RGC to 0.754% for rods, which provides no evidence for a cell-type-specific
filtering artifact.

Donor depth imbalance is greatest for AII amacrine (CV 0.511), astrocyte (0.472),
microglia (0.406), RGC (0.387), and rods (0.346). All retain usable depth, but
these values reinforce the requirement for sample-specific library-size
normalization and donor-blocked contrasts.

The ratio of unique fragment rows to fifth-column read support shows a highly
consistent donor pattern across every cell type: mean ratios are 0.669 for LGS1,
0.680 for LGS2, 0.725 for LGS3, and 0.613 for LVG1. This is a donor/library
property rather than a cell-type-specific anomaly. Differential testing must use
unique fragment rows, donor blocking, and estimated normalization factors; the
read-support column must not be expanded into duplicate inferential counts.

Peak calling should be performed for each donor and for each pooled cell type.
A pooled cell-type peak should be considered reproducible only if it overlaps a
peak in at least two of the four donor pseudobulks. Standardize reproducible
summit-centered peaks to a common width, remove blacklist and noncanonical
regions, and resolve overlaps deterministically across all 13 cell types. Count
all 52 pseudobulks over the resulting common universe and review FRiP, library
sizes, peak counts, sample correlations, donor effects, and cell-type separation
before differential accessibility.
