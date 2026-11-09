# Top 30 natural MG parents: visual peak atlas

The portable seven-page PDF is at
`output/pdf/MG_top30_parent_peak_atlas.pdf` relative to the project root.

## Read the figures

- **Page 1:** 30 regions by cell type. The left heatmap shows observed mean
  browser signal over the exact 500bp parent. One linear color scale is shared
  across every region and cell type; there is no row normalization. The middle
  heatmap shows the original donor-blocked pairwise log2FCs. The right column
  gives each region's largest pairwise FDR.
- **Pages 2-7:** five regions per page, each compared across all 13 cell types.
  Read down a column. Every track within that column has exactly the same y-axis
  range, printed above the MG track. There is no independent autoscaling by cell
  type. Different loci can use different ranges. Black is MG, orange astrocytes,
  blue microglia, and gray the neuronal types.
- The gray vertical shading is the exact 500bp natural parent. The complete
  window spans 2.5kb. Neighboring peaks outside the shading are not the parent
  insert. Coordinates printed above tracks use IGV conventions: hg38,
  1-based inclusive. Parent numbers refer to original candidate ranks.
- Parenthetical gene labels identify the nearest annotated gene. They are
  descriptive and do not establish a regulatory target gene.

The matching high-resolution PNG files can also be opened individually.

## Selection rule

Start from the existing 1,000 parents in
`reports/phase3/candidate_ranking/top_1000_candidates.tsv.gz`. For every region,
take the **largest original FDR among its 12 MG-versus-cell-type comparisons**.
Sort this value ascending, using original candidate rank to break ties, and take
the first 30. This emphasizes evidence against the hardest comparison rather
than allowing a strong MG-versus-neuron result to conceal weak glial evidence.
The maximum is a selection metric, not a newly calculated combined FDR.

The rule was fixed before inspecting track shapes. This is not the original
composite top-30 candidate ranking, nor a selection of the largest fold changes.
Existing FDRs retain their original adjustment over all 147,359 tested regions
within each contrast. No differential model was rerun and no multiple-testing
correction was recalculated on this subset.

These 30 have largest pairwise FDRs from 1.19e-9 to 2.50e-7. The smallest fitted
pairwise log2FC across all 30 and 12 comparisons is 1.576, approximately 2.98-fold.
Mean observed browser MG signal divided by the strongest off-target ranges from
2.76 to 10.27, with median 4.15. All 30 retain FDR < 0.05 in all eight glial
leave-one-donor-out refits.

## Interpretation and limitations

The figures support stronger observed MG accessibility. Astrocyte and microglial
peaks are often still visible. These are selected best statistical examples;
they are not representative of all 1,000 parents or independent validation.
Natural accessibility and differential significance do not establish enhancer
activity in a reporter or specificity of designed variants.

Browser tracks use observed fragment endpoints, pooled global insertion scaling,
and 100bp mean bins. DA uses fragment overlap counts, sample-specific TMM
normalization, and biological donor blocking. The two fold measurements are
therefore not expected to be identical. Low off-target glial depth limits the
interpretation of absent tracks. All four donors contributed to the analysis.

## Reproduce and audit

```bash
MPLCONFIGDIR=/tmp/glia-matplotlib .venv/bin/python scripts/plot_top30_parent_peak_atlas.py
```

- `top30_selected_parents.tsv`: identities, coordinates, all original DA fields,
  display order, selection statistic, and browser ratios.
- `top30_browser_signal.tsv`: exact interval means for all 13 types.
- `displayed_track_values.npz`: the plotted 2,500 values for each of 30 regions
  and 13 cell types, with group order and relative coordinates.
- `selection_and_provenance.json`: source hashes, method, parameters and summary.
- `validation.json`: numerical and PDF validation results.

All 13 pooled track hashes were checked against the delivered package manifest.
All 390 parent browser means were checked against the previous whole-library
audit. Every displayed track uses the same y-axis as the other cell types at its
locus, and every plotted value is finite and nonnegative. All seven final PDF
pages were rendered and visually inspected.
