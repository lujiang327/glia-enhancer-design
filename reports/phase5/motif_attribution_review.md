# Motif/parent attribution review

Job **63285123** completed successfully (exit 0:0) in **61 seconds**, with maximum RSS 5,326,092KB. All four downloaded output checksums passed. The review independently verified preservation of all original library and motif comparison evidence, attribution ranges, orientation sign flags and single-base edit direction. All **45 tests** pass.

The library still contains **1,000 natural parents, 938 gain variants and 946 loss controls**. Attribution context was added to **59,828 edit-overlapping motif comparisons**. The complete library is `motif_attribution/ranked_library_with_parent_attribution.tsv.gz`; candidate-specific motif evidence is in `edit_motif_attribution_review.tsv.gz`. Original DA statistics, natural sequences, variant sequences, direct predictions, repeat context and reference mappability remain preserved.

| Measure | Gain variants | Loss controls |
|---|---:|---:|
| Edits | 938 | 946 |
| Count hypothetical edit direction agrees in both orientations | 938 | 946 |
| Median edited-base absolute parent count-attribution percentile | 90.6% | 99.8% |
| Median direct mean logcount change | +0.584 | −1.255 |
| Median conservative effect magnitude across orientations | 0.478 | 1.089 |
| At least one non-enzyme match window with positive count attribution in both orientations and a top-10% parent base | 670 | 940 |
| A lost match in that positive parent context | 459 | 919 |
| A lost match in negative parent context with a top-10% parent base | 231 | 3 |
| Hypothetical attribution delta versus direct model effect, within-class Spearman | 0.297 | 0.665 |

Loss controls mostly perturb bases with very high parent count attribution and remove a motif match in a positively attributed interval. This is consistent with disrupting learned positive sequence features. The gains have weaker attribution-effect correlation and mixed parent motif context; some lose matches in positive windows while increasing predicted counts. Therefore, a gained/lost motif label or the sign of its whole window is insufficient to explain an edit. A window can span both positive and negative bases, and overlapping motif families share the same positions.

The all-edit direction agreement is a **consistency check, not independent validation**: these same parent attribution maps were used to nominate edits, and the retained edits were screened for direct model effects. The within-class correlations likewise describe selected edits and are not an unbiased benchmark of all possible substitutions. Do not compare hypothetical attribution values numerically as calibrated direct logcount changes.

There are **360 enzyme-motif comparisons** (explicit DNASE_/TN5_ database entries). They are retained and flagged, but excluded from the descriptive non-enzyme context counts above. Other matches remain database motif identities; no TF expression, binding, regulatory role or causal family identity is assigned automatically. Example cases follow original DA rank, selecting five high-attribution windows per edit from the first ten gains and first ten controls; enzyme identities are flagged. These are review examples, not a newly selected synthesis set.

All site attribution describes the **natural parent**, including intervals where a motif appears only after mutation. Saved float16 maps introduce quantization. Profile attribution concerns profile shape, and its sign is not a total-accessibility criterion. The model scores an insert in its native 2,114bp context; isolated reporter behavior can differ. The current figure shows selected edited-base percentiles and direct versus hypothetical effects and is stored as `motif_attribution/motif_attribution_review.png`.

**Recommendation:** proceed to experimental library planning and manual review of representative parent/gain/control sets. No new failure in attribution consistency was found, and no additional hard exclusion is introduced here. Continue to label these as model-guided accessibility variants, not validated enhancers or demonstrated specificity improvements. Off-target variant prediction, causal motif validation and reporter performance remain unresolved. Final insert/adaptor/barcode/cloning constraints and the number of experimental constructs must be specified before a synthesis order.

Reproduce verification and summaries from repository root:

```bash
MPLCONFIGDIR=/tmp/glia-matplotlib .venv/bin/python scripts/review_phase5_motif_attribution.py \
  --results-dir reports/phase5/motif_attribution
```

The source checksum manifest remains unchanged and covers downloaded outputs only. `review_output_checksums.sha256` covers the derived review files separately. Source job, code revision, container and attribution hashes are preserved in the downloaded `runtime.txt` and `motif_attribution_summary.json`.
