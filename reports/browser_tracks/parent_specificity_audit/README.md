# Review of natural-parent MG enrichment

Reviewed October 20, 2026. This analysis checks existing results and the delivered
browser package. It does not refit differential tests, train models, or change
the selected library.

## Findings

- All 1,000 natural parents have positive donor-blocked MG effects and original
  FDR < 0.05 against each of the 12 off-target cell types. The largest pairwise
  FDR is 0.00840. These FDRs were adjusted across the original 147,359 tested
  regions separately for each contrast, not recomputed within the selected set.
- MG-versus-astrocyte median log2FC is 2.447 (about 5.45-fold fitted enrichment),
  with median FDR 8.50e-6. MG-versus-microglia median log2FC is 2.035 (about
  4.10-fold), with median FDR 8.60e-5. Related glia are the closer comparisons.
- Direct reading of the delivered 100bp-binned BigWigs gives higher mean MG
  signal than every off-target for all 1,000 exact 500bp parent intervals.
  MG divided by the highest of 12 off-target means ranges from 2.258 to 11.854,
  with median 3.837. Only 426 parents have a ratio at least four. Signal in
  off-target cells is often clearly present.
- 978 parents have MG signal above both astrocytes and microglia in each of
  the four donor browser tracks. The other 22 fail this descriptive comparison
  for at least one donor. Per-donor passing counts are LGS1 981, LGS2 997,
  LGS3 1,000 and LVG1 1,000. This comparison is not a significance test.
- 872 parents retain FDR < 0.05 in all eight glial leave-one-donor-out refits.
  The other 128 lose significance in at least one refit. All selected parents
  retain positive effects in those refits. Lower power and donor influence
  cannot be fully separated with four donors.
- All 25 delivered BigWig SHA256 hashes match the package manifest. All 1,000
  candidate BED intervals and identifiers match the DA table. This excludes a
  mismatch among the inspected package files, not all possible upstream errors.

## How to review the files

The Excel workbook is
`outputs/parent_specificity_review_20261020/MG_1000_parent_specificity.xlsx`
relative to the repository root. Start with **Review**, then **Parents**. The
**log2FC** and **FDR** tabs show one column per comparison, retaining all 1,000
rows. **Browser signal** has observed mean insertion density for all 13 cell
types. Scientific notation prevents small FDRs from appearing as rounded zeros.
In Parents, filter worst glial LODO FDR >= 0.05 or browser donor count < 4 to
inspect sensitivity exceptions. Rank and peak ID link every tab.

`specificity_audit.png` summarizes effects, original FDRs, and observed browser
ratios over the entire selected set. `shared_scale_locus_examples.png` shows
fixed ranks 1, 500 and 1,000 plus rank 972, the weakest pooled browser ratio.
Each locus uses the same y-axis across its three cell-type tracks. Different
loci can use different scales. The exact 500bp parent is shaded.

Rank 1 is `chr13:27504893-27505392` in IGV coordinates. Rank 972 is
`chr2:121821205-121821704`. View roughly 1kb of flanking sequence on each side.
In IGV select hg38, load the candidate BED and relevant pooled BigWigs,
multi-select the quantitative tracks, and choose **Group Autoscale**. Independent
autoscaling can make lower-amplitude off-target peaks appear as tall as MG.
Zooming far out can also obscure the candidate interval among nearby peaks.
See the [official IGV quantitative-track instructions](https://igv.org/doc/desktop/UserGuide/tracks/quantitative_data/).

We have not inspected Thanh's exact IGV locus or settings. Scaling is a possible
explanation for the visual impression, not an established diagnosis of his view.

## What these numbers establish

The appropriate description is **MG-enriched accessible regions**. There are
shared off-target peaks, particularly in related glia. Statistical significance
does not require an off-target signal to be zero. The regions were selected for
specificity using these data, so favorable enrichment is expected and is not
independent validation.

Browser tracks count fragment endpoints, use pooled whole-track insertion
normalization and weight donors by depth. DA counts fragments overlapping each
consensus interval, uses sample-specific TMM normalization, and blocks on
biological donor. Thus a pooled browser ratio need not equal 2 raised to the
fitted DA log2FC. Both measurements are compositional, not absolute accessibility
per cell. Off-target glial depth is limited, and a missing peak is weak evidence
of absence. The 100bp smoothing also limits interpretation of fine structure.

None of these measurements establishes reporter activity, exclusive MG activity,
or preserved specificity of designed variants. No off-target sequence models
were evaluated by this review. Before interpreting a region as unsuitable, inspect
its exact parent interval, effect size, FDR, and donor behavior together. Before
ordering a large specificity-focused library, the previously identified variant
and experimental evidence gaps still apply.

## Reproduction and numerical exports

```bash
MPLCONFIGDIR=/tmp/glia-matplotlib .venv/bin/python scripts/audit_parent_specificity.py
```

Input DA table: `reports/phase3/candidate_ranking/top_1000_candidates.tsv.gz`.
Input browser folder: `data/intermediate/browser_tracks/retina_hg38_review_100bp`.
`source_checksums.json` records source identities.

TSV exports preserve the source statistical precision and give flat tables for
other tools: `parents_review.tsv`, `log2FC.tsv`, `FDR.tsv`, `browser_signal.tsv`,
`donor_browser_signal.tsv`, and `contrast_summary.tsv`. Coordinate start/end
columns use BED conventions (0-based, end exclusive). The IGV locus column is
1-based inclusive.

`scripts/create_parent_specificity_workbook.mjs` builds the workbook using the
bundled artifact-tool runtime and `workbook_inputs.json`. Copy the builder to a
temporary directory with a `node_modules` link to the runtime dependencies, then
run it with the audit directory and workbook output directory as arguments.
Workbook summaries link to the source metric sheets with formulas. XML export
checks confirm all 13,000 FDR values and 13,000 log2FC values against the source.
