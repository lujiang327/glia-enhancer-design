# GSE196235 donor-by-cell-type audit review

**Decision: GO for donor-by-cell-type pseudobulk and peak preparation, with
caveats. Differential accessibility remains blocked until those outputs pass
depth, library-complexity, peak-universe, and donor-consistency QC.**

All 13 published cell-type fragment pools contain unambiguously assigned cells
from all four biological donors. The retained data comprise 48,422 cells,
781,336,301 canonical chr1-22/X fragment rows, and 1,172,284,901 canonical
fifth-column read-support counts. Depending on cell type, 92.67-95.67% of bare
barcodes, 91.40-94.89% of fragment rows, and 91.54-94.85% of read support are
retained by the unambiguous donor mapping. This loss is consistent across cell
types and does not indicate a cell-type-specific mapping failure.

Muller glia contribute 3,434 cells and 86,606,259 canonical fragment rows across
LGS1, LGS2, LGS3, and LVG1. All major neuronal off-targets have four-donor
coverage and substantial depth. AII amacrine also has at least 59 cells per
donor and 10.6 million canonical fragment rows in total.

Four donor-by-cell-type strata contain fewer than 50 cells:

| Cell type | Donor | Cells | Canonical fragment rows | Canonical read support |
|---|---:|---:|---:|---:|
| Astrocyte | LGS1 | 12 | 437,180 | 660,271 |
| Astrocyte | LGS2 | 46 | 1,730,821 | 2,545,694 |
| Microglia | LVG1 | 27 | 666,219 | 1,077,671 |
| RGC | LGS1 | 42 | 1,770,885 | 2,629,060 |

These strata remain usable for counting over a common peak universe, but their
cell numbers make cell-type-specific peak recovery and pairwise inference less
stable. Retain them in the primary dataset, flag the corresponding pairwise
contrasts, require consistent effect direction across donors, and run a
sensitivity analysis that omits the low-cell donor for the affected comparison.

Rod photoreceptors contribute 30,062 cells and 380,954,728 canonical fragment
rows, far more than any other cell type. Do not form a raw pooled "rest" library,
because rods would dominate it. The primary Muller-versus-rest contrast must use
an explicit equal-weight cell-type contrast or another prespecified method that
prevents abundance from determining off-target weight. Retain the planned paired
Muller-versus-each-cell-type contrasts.

Next, construct the 52 donor-by-cell-type pseudobulks while retaining both unique
fragment-row depth and fifth-column read support. Call peaks from pooled
cell-type data, create a blacklist-filtered canonical consensus union across the
target and off-target cell types, and count every donor pseudobulk over that same
universe. Review pseudobulk depth, FRiP, library complexity, peak overlap, sample
correlations, and donor/cell-type separation before enabling differential testing.
