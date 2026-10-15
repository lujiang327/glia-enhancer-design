# MG browser-review BigWigs

The browser package is generated under ignored `data/intermediate/browser_tracks/MG_hg38_100bp/`. Large binary signal files remain outside Git. A matching tar archive is prepared alongside this directory for manual sharing; nothing is sent or published automatically.

## Contents

- `MG_pooled.IPM_per_base.100bp.bw`: observed pooled Müller-glia ATAC insertion density.
- `MG_LGS1`, `MG_LGS2`, `MG_LGS3`, `MG_LVG1` `.IPM_per_base.100bp.bw`: the same signal separately for each biological donor.
- `MG_top1000_parents.hg38.bed`: all ranked natural 500bp candidate regions; feature names contain rank and peak ID.
- `candidate_index.tsv`: ranks, coordinates and genomic/nearest-gene annotations.
- `first_locus.txt`: a browser-ready interval around the highest-ranked candidate.
- `browser_track_manifest.json`, `hg38.canonical.chrom.sizes` and `checksums.sha256`: signal definition, source hashes, chromosome lengths and output hashes.

All files use hg38, chromosomes 1–22/X. The large raw single-base insertion BigWigs remain unchanged at `data/intermediate/training/muller.insertions.bw` and the four donor `.insertions.bw` paths. Use those when single-base structure is needed; the overview files cannot resolve individual insertion footprints.

## Signal definition

For each fixed 100bp bin (shortened at a chromosome end):

`signal = bin insertion count / actual bin width × 1,000,000 / retained track insertion count`.

Units are **insertions per million retained insertions per base**, averaged over 100bp. Missing source intervals represent zero insertions; output zero bins are explicit. Each track's signal integrated over the genome is approximately one million, with float32 encoding/BigWig-header rounding tolerance. Each fragment contributes its two previously validated insertion endpoints once. The fifth-column read-support count is not used as a multiplicity weight. No new Tn5 endpoint shift is applied.

The pooled raw track has 171,959,170 insertions from 85,979,585 retained fragments. Every donor is normalized by its own insertion total; the pooled track derives from pooled counts, so it is depth-weighted rather than an equal-donor mean. Global library normalization does not provide absolute accessibility per cell and does not remove all composition or batch effects.

## Open in IGV

1. Select the **Human hg38/GRCh38** reference genome.
2. Load the five `.bw` files and `MG_top1000_parents.hg38.bed` as local data files.
3. Paste the interval from `first_locus.txt` into the locus search box, or use coordinates from `candidate_index.tsv`.
4. Multi-select the numeric tracks, right-click and choose **Group Autoscale** to compare their normalized signal using a common vertical scale. Independent autoscaling can make unequal peaks appear similar.
5. Keep the parent BED track visible to distinguish the exact 500bp candidate from surrounding accessible chromatin.

See the official [IGV quantitative-track instructions](https://igv.org/doc/desktop/UserGuide/tracks/quantitative_data/) and [reference-genome instructions](https://igv.org/doc/desktop/UserGuide/reference_genome/).

## Interpretation and missing comparisons

These are **observed natural-locus ATAC tracks**, not model predictions, RNA expression, enhancer assays or synthetic-variant measurements. MG donor agreement is useful evidence about the natural parent locus. Off-target tracks are not included in this first package. To review specificity, we also need matched normalized astrocyte/microglia and relevant neuronal signal tracks; those can be built from the existing donor pseudobulks without training off-target sequence models.

The reusable Great Lakes SSH session was unavailable during this export. Local MG inputs were sufficient to finish this package; remote off-target track inventory/preparation requires a renewed session. Do not interpret the MG-only view as evidence that other cell types are inactive.

## Reproduce

Run from the repository root, using a new output directory:

```bash
.venv/bin/python scripts/prepare_browser_tracks.py \
  --output-dir data/intermediate/browser_tracks/MG_hg38_100bp
```

The exporter verifies raw insertion totals, conserves integrated counts during binning, checks normalized signal integrals and records hashes. Unit tests cover boundary-spanning intervals, chunk-boundary clipping, zero intervals, partial terminal bins and a real BigWig round trip. It refuses to overwrite an existing output directory.
