# Matched retinal cell-type browser tracks

The renewed Great Lakes session confirmed the existing four-donor pseudobulks are available. No existing off-target BigWigs were found in the inspected Phase 3 scratch tree. This workflow creates observed ATAC tracks for comparison with the MG-only package; it does not train or evaluate off-target sequence models.

## Outputs and resources

The CPU job array covers MG, astrocytes, microglia and all ten neuronal off-target cell types. Each cell type gets a pooled hg38 100bp BigWig. MG, astrocytes and microglia additionally get four separate donor BigWigs. The full result has 25 tracks: 13 pooled and 12 donor-resolved. A compact initial comparison can use MG, astrocyte and microglia pooled tracks, adding neuronal tracks and donor lanes as needed.

The array requests one CPU and 4GB RAM per task, with a six-hour maximum and at most three concurrent tasks. The runtime limit is a ceiling, not an estimated duration. No GPU is requested. Outputs go to `browser_tracks_100bp/<cell_type>/` under the existing private Phase 3 scratch root. The account and scratch path are supplied manually, not stored in tracked settings.

## Input and signal contracts

`config/browser_celltype_tracks.json` pins all 52 donor pseudobulk hashes, expected fragment totals and cell counts from the reviewed Phase 3 QC, as well as chromosome sizes and the original raw MG track hash. Every input checksum is verified before use. Fragment rows are processed once, never expanded by their read-support column. Donor identity and low-cell strata remain visible in the summary.

Insertions are placed at the same `start` and `end` endpoint coordinates used for the existing MG track, with no new shift. Coordinates outside the represented reference fail rather than being silently dropped. The MG array task rebuilds integer 100bp counts and requires exact genome-wide agreement against the independently prepared original MG single-base BigWig. This is a consistency check of the export convention; it does not independently estimate Tn5 bias for every cell type. Require this MG check to pass before combining the exported tracks for review.

All cell types use the MG browser-package units: `bin insertions / actual bin width × 1e6 / retained track insertions`. Each donor uses its own insertion total, and the pooled track uses pooled counts. Output bins explicitly include zeros. Each normalized BigWig must integrate to approximately one million; input and binned insertion counts must equal exactly twice the retained fragment count.

Files remain ignored large intermediates. Compact manifests, QC and previews can be downloaded into reports after completion, then a selected sharing package can combine observed cell-type tracks with the existing 1,000-region BED file. The original single-base training tracks are not overwritten. Existing output directories are never overwritten; failed jobs retain staging inputs for diagnosis.

## Submit after the user pushes

From the Great Lakes repository root, after pulling the committed workflow and exporting the private scratch root:

```bash
sbatch --account="$SLURM_ACCOUNT" hpc/slurm/build_celltype_browser_tracks.sbatch
```

Monitor the returned array job ID with `squeue` and `sacct`; each task writes `logs/slurm/retina_browser-<array_job>_<task>.out/.err`.

## Interpretation

Use the same hg38 assembly, candidate intervals and shared vertical scale in IGV. Tracks allow visual inspection of natural-locus signal and donor variation. Lower-depth astrocyte/microglia tracks may be noisy: a small or absent peak is not evidence that accessibility is truly absent. Global library normalization is compositional, not absolute accessibility per cell. Pairwise donor-blocked DA statistics remain the inferential evidence; browser views supplement them. None of these observed tracks measures synthetic-variant activity or establishes preserved mutant specificity.

Tests cover four-donor pooling, endpoint positions, the shortened terminal bin, normalized signal, a real raw-MG BigWig reference comparison, count/coordinate mismatches and refusal to overwrite outputs.
