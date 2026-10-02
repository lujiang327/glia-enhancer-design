# Phase 4 natural-locus annotation protocol

## Scope

Annotate all 1,000 original 500-bp natural parents and the 1,884 provisional
best passing single-base edits. Preserve the original rankings, sequence
predictions, and observed specificity evidence. This stage adds genomic
context and makes no automatic repeat or mappability exclusions.

## Reference sources

Use the hg38 `rmsk` table from UCSC RepeatMasker and retain the accompanying
SQL schema. The parser verifies its 17-column schema, uses genomic chromosome,
start, and end fields, and retains repeat name, class, and family. It counts
the union of clipped repeat intervals rather than summing overlapping records.
The schema and dataset are available from the [UCSC hg38 annotation database](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/).

Use the Umap **multi-read** bigWigs for unconverted hg38 at read lengths 50
and 100 bp, available from the [UCSC Hoffman mappability archive](https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/).
These lengths provide two descriptive reference-locus sensitivity views.
They do not reproduce paired-end ATAC alignment with the actual read lengths,
aligner, mismatch tolerance, or read-quality filters.

The [Umap source documentation](https://bismap.hoffmanlab.org/) defines
multi-read mappability and states that zero values are omitted from its
bigWig files. Therefore this workflow converts omitted positions to zero
and divides by the full 500-bp parent length. It reports the fraction of
omitted positions separately. Chromosome names and lengths must match the
existing hg38 reference FASTA index; missing chromosomes cause failure.

## Provenance and storage

The CPU job downloads the public reference assets into ignored
`data/raw/phase4_annotation/`. The two mappability files together occupy
approximately 2 GB; RepeatMasker and its schema add further storage.
The first download records the source and resolved URLs, retrieval time,
HTTP Last-Modified/ETag where available, byte count, and observed SHA256.
These are observed snapshot hashes rather than upstream published checksums.
Subsequent runs verify the saved hashes and refuse changed or unmanifested
assets. Full downloads are preserved for later reproducibility.

The job runs on the `standard` CPU partition with 2 CPUs, 8 GB RAM, and a
4-hour limit. It verifies container dependencies before downloading, checks
the source snapshot before annotation, and saves outputs under the private
scratch root's `phase4_locus_annotation/`. Raw tracks remain in the project
raw-data directory and are excluded from Git. Account and private directory
configuration are supplied at submission rather than committed.

## Outputs and interpretation

- `annotated_natural_parents.tsv.gz`: all original parent data plus union
  repeat overlap, intersecting repeat names/classes/families, and 50/100-bp
  Umap mean, minimum, zero fraction, fraction at least 0.9, and omitted fraction.
- `annotated_best_single_edits.tsv.gz`: the provisional edits with their
  parent annotations, edited-base repeat context, and reference-locus Umap
  values at the edited position.
- `parent_repeat_hits.tsv.gz`: repeat coordinates and labels for auditing
  parent overlaps.
- `annotation_summary.json`, `annotation_provenance.json`, `runtime.txt`,
  and `output_checksums.sha256`: summary, asset/input fingerprints, and logs.

An edited-base Umap value describes the **natural reference sequence** at that
position. This workflow does not recompute genomic uniqueness for the variant.
Neither a repeat overlap nor a low Umap score alone proves a region is unusable.
Repeat-derived regulatory activity is possible, and the mapping score is a
technical context metric. Review these annotations with accessibility,
donor consistency, model behavior, and assay design before selecting oligos.

Off-target accessibility of designed sequences remains unavailable; retain
the distinction between observed parent specificity and predicted variant
specificity. After annotation review, assemble the parent/variant library
table with the differential statistics, important motifs, prediction changes,
and donor-reproducibility evidence.
