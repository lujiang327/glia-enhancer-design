# Activity-first top 30 with MG-specificity safeguards

Prepared October 26, 2026. This is a review-only re-ranking of the fixed 1,000
natural parents. The original candidate library and significance atlas are
preserved. No model was trained, differential test rerun, or sequence changed.

## Files to share

- `top30_activity_specificity.csv`: the compact ranked list, readable in Excel.
  It includes original rank, peak ID, hg38 coordinates, nearest gene, observed
  MG signal, equal-donor TMM MG CPM, specificity, and donor/statistical metrics.
- `output/pdf/MG_top30_activity_specificity_atlas.pdf` from the repository root:
  seven-page atlas showing all 30 regions against all 12 off-target cell types.
- `activity_vs_significance.png`: direct comparison with the previous atlas.

The peak atlas now has **one global y-axis range across all 30 regions and all
13 cell types**. Shading marks the exact 500bp parent; each window spans 2.5kb.
The numeric signal is globally normalized observed ATAC insertion density.
No cell type or region is independently rescaled. Gene names are nearest-gene
annotations, not established regulatory targets.

## How the original atlas was ranked

For region i, define W_i as the maximum of its original donor-blocked
MG-versus-cell-type FDRs across the 12 off-targets. Sort W_i ascending and select
30; break ties by the original candidate rank. Smaller W_i means the weakest
comparison is still statistically well supported. It is a ranking statistic,
not a new combined FDR. FDRs retain the original multiple-testing adjustment.

This selects consistent differential evidence. It does not maximize measured
MG accessibility, peak height, ChromBPNet prediction, or reporter activity.
The original 1,000-parent composite ranking, the previous atlas ranking, and
this new activity-first ranking are three distinct orders.

## New activity-first ranking

For each exact 500bp parent, read the delivered 100bp-binned BigWigs at their
full stored resolution, including partial bins:

- A_i = mean MG browser signal over the parent.
- B_i = largest mean signal among the 12 off-target cell types over that parent.
- S_i = A_i / B_i. No pseudocount is needed: B_i is positive for all 1,000.

Eligibility requires:

1. S_i >= 3, so observed mean MG signal is at least threefold above every
   off-target's mean signal in the interval.
2. FDR < 0.05 and fitted log2FC >= 1.5 in every pairwise comparison. This is
   approximately 2.83-fold minimum fitted enrichment; it is a separate
   donor-blocked/TMM measure from the pooled browser ratio.
3. FDR < 0.05 in all eight astrocyte/microglia leave-one-donor-out refits.
4. Positive MG-versus-rest effects in all four donors and MG browser signal
   above both other glia in every donor's tracks.

787 of the 1,000 parents pass. Within this set, sort A_i descending, then S_i
descending, then original candidate rank ascending. Take 30. These are explicit,
adjustable review choices, not universal biological or assay-validation cutoffs.
The ranking targets interval-wide accessibility, rather than just the tallest
100bp bin. All component values and the full 1,000-region eligibility table are
retained for auditing.

## What changed

| Measurement | Previous significance top 30 | New activity-first top 30 |
|---|---:|---:|
| Median MG browser density, IPM/base | 0.03176 | 0.03612 |
| Median equal-donor TMM MG CPM | 42.96 | 46.51 |
| Minimum MG / strongest off-target browser ratio | 2.76 | 3.03 |
| Median MG / strongest off-target browser ratio | 4.15 | 3.55 |

The new list increases median measured MG browser signal by **13.7%**, but its
median specificity ratio decreases. Fourteen parents overlap. Every new parent
has all pairwise FDRs below 9.13e-6 and passes all eight glial donor-omission FDR
checks. Those same checks helped select the list; they are not independent
validation.

The improvement is modest because the previous list already contained several
of the highest-accessibility parents. Even taking the 30 largest MG browser
signals among all 1,000 without any additional specificity or donor filters
raises the top-30 median by only **20.4%**. Thus an order-of-magnitude increase
cannot be obtained by re-ranking these fixed 1,000 under this signal measure.
To investigate a materially different activity range, inspect the broader
previously scored candidate universe rather than merely changing the top-30
ordering. That expansion has not been performed here.

## Specificity versus accessibility sensitivity

Holding the other new safeguards fixed and changing only the observed browser
ratio floor gives:

| Required browser ratio | Eligible parents | Top-30 median MG density | Top-30 median ratio |
|---|---:|---:|---:|
| >=2.5 | 833 | 0.03825 | 3.47 |
| >=3.0 (main list) | 787 | 0.03612 | 3.55 |
| >=3.5 | 579 | 0.03157 | 3.89 |
| >=4.0 | 354 | 0.02699 | 4.48 |

A stricter fourfold requirement yields a more specific but lower-accessibility
top 30. `threshold_sensitivity.csv` records the results so this tradeoff remains
visible. None of these ratios is evidence of MG-exclusive reporter activity.

## Interpretation and limitations

Browser signal is pooled global insertion scaling, in insertions per million
insertions per base (IPM/base), averaged over 100bp bins. Small numeric values
are a consequence of the units and do not by themselves establish a weak peak.
These are compositional, depth-normalized measurements, not absolute
accessibility per cell. The equal-donor TMM CPM column is a useful companion
measure; it is not interchangeable with the pooled density.

The ranking is restricted to the existing parents and uses the same observed
data that selected them. It is not independent validation, a ChromBPNet
sequence score, or a prediction of vector expression. The broader annotation,
rod-RNA, off-target-model and experimental-validation caveats remain. Sequence
and assay suitability must still be reviewed before ordering.

## Reproduce

Parameters are saved in `config/activity_specificity_top30_review.json`.

```bash
MPLCONFIGDIR=/tmp/glia-matplotlib .venv/bin/python scripts/plot_top30_parent_peak_atlas.py \
  --ranking activity_specificity \
  --output-dir reports/browser_tracks/top30_activity_specificity \
  --pdf output/pdf/MG_top30_activity_specificity_atlas.pdf
```

Source DA table: `reports/phase3/candidate_ranking/top_1000_candidates.tsv.gz`.
Source browser measurements: `reports/browser_tracks/parent_specificity_audit/`.
The atlas reopens all 13 pooled BigWigs, checks their SHA256 hashes, and compares
the selected interval means against the whole-library audit. Full tables,
plotted values, original ranks, threshold sensitivity, parameters and hashes are
retained. All seven PDF pages were rendered and visually inspected.
