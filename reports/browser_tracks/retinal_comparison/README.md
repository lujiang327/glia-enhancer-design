# Retinal cell-type browser comparison

The source array is **64054545**, with 13 cell-type tasks. Each task produced its observed browser tracks and exited 0:0. The original Slurm accounting is saved in `../retinal_array_accounting.psv`. The separate short recovery allocation **64072222** was used only to access/stage files when login-node scratch reads stalled; it did not regenerate the scientific tracks.

**Verified result:** all 25 downloaded BigWigs passed checksums, donor/input
provenance, chromosome and normalization checks. The MG reference comparison
passed exactly across the genome. All **54 tests** pass. The predefined six-locus
preview shows stronger MG signal within the parent regions than the two glial
off-targets, but these selected natural regions are not independent validation of
mutant activity or specificity.

The assembled sharing directory is kept outside Git at `data/intermediate/browser_tracks/retina_hg38_review_100bp/`. Its tar archive is prepared alongside it for manual sharing. No files are sent or published automatically.

- **Compact first look (~418 MB):** `retina_hg38_review_100bp_glial_core.tar`,
  containing the three pooled glial tracks, candidate BED/index and preview.
- **Full set (~3.43 GB):** `retina_hg38_review_100bp.tar`, containing all 25 tracks.

Both archives have their own file-checksum lists; `sharing_archives.json` records
the archive SHA256 values. Archives omit machine-specific owner metadata. The
compact subset does not include neuronal or donor tracks; review those from the
full package when assessing broader specificity and donor consistency.

## Files to inspect first

Start with these three pooled tracks plus `MG_top1000_parents.hg38.bed`:

- `MG_pooled.IPM_per_base.100bp.bw`
- `Astrocyte_pooled.IPM_per_base.100bp.bw`
- `Microglia_pooled.IPM_per_base.100bp.bw`

All ten neuronal off-target pooled tracks are also present. The full set includes 25 BigWigs: 13 pooled cell types and four donor tracks each for MG, astrocytes and microglia. `track_index.tsv` identifies every file, its cell type, donor, insertion total, byte size and SHA256. `candidate_index.tsv` contains all 1,000 ranked natural regions and their annotations.

Use **hg38/GRCh38** in IGV, load the `.bw` files and candidate BED, and select **Group Autoscale** across tracks you wish to compare. Independent autoscaling can exaggerate similarity between peaks with different heights. Keep the exact 500bp parent interval visible; an adjacent peak is not necessarily within the tested insert. The pooled tracks are depth-weighted over donors rather than equal-donor means.

`MG_vs_glial_offtargets_preview.png` shows prespecified ranks 1, 200, 400, 600, 800 and 1000 in all three pooled glial tracks, with the same vertical scale across every panel. These ranks were not selected for visual peak quality. Donor tracks let reviewers inspect whether a pooled pattern is carried by one donor. Astrocyte donors LGS1/LGS2 and microglia donor LVG1 have especially few cells, so absence of a visible peak is weak evidence of biological absence.

## Signal and interpretation

All tracks use observed native-locus fragment endpoints and the same 100bp mean signal: `bin insertions / actual bin width × 1e6 / whole-track retained insertion count`. Zero bins are explicit. Each input fragment row contributes two endpoints once; read-support counts do not inflate molecule counts. The rebuilt MG track must match the original validated raw MG BigWig exactly after integer binning across the genome.

The assembler verifies original output checksums, group/config identity, all donor hashes/cell counts/fragment totals, track inventories, chromosome lengths, normalized integrals and the MG reference result before creating the sharing folder. Original source summaries and compact review manifests are kept with this report. Large binaries remain outside Git.

These are descriptive browser views, not formal differential tests or model outputs. Global library normalization differs from the sample-specific edgeR normalization used in donor-blocked DA, so visual height ratios need not equal the DA log-fold changes. Normalization is compositional and does not measure absolute accessibility per cell. Lower-depth cell types have noisier profiles. Use the original donor-aware statistics alongside the browser view.

**No off-target sequence model was trained.** These tracks show natural-locus accessibility only; they do not predict whether a designed gain mutation also increases accessibility in another cell type and do not establish reporter activity.

## Reproduce after downloading completed source outputs

```bash
MPLCONFIGDIR=/tmp/glia-matplotlib .venv/bin/python scripts/assemble_retinal_browser_package.py \
  --tracks-root data/intermediate/browser_tracks/retinal_verified_outputs \
  --output-dir data/intermediate/browser_tracks/retina_hg38_review_100bp \
  --report-dir reports/browser_tracks/retinal_comparison
```

The output directory must be new. The assembler uses hard links where possible to avoid duplicating local BigWig bytes; tar archives contain the actual file content and remain portable when extracted. Original source checksums and a separate final sharing-package checksum list are retained.

Official [IGV quantitative-track instructions](https://igv.org/doc/desktop/UserGuide/tracks/quantitative_data/) describe Group Autoscale.
