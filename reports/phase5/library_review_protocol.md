# Draft parent/variant library and candidate motif scan

The draft retains all 1,000 ranked natural 500bp parents and the best passing single-base edits: 938 gain variants and 946 loss controls, 2,884 sequences total. Rows are ordered by original differential-accessibility candidate rank, then parent, gain, control. This is a review library, not a synthesis order or a new ranking based on variant specificity.

`ranked_library_draft.tsv.gz` combines donor and differential-accessibility statistics, genomic annotations, natural-locus repeat/mappability context, parent attribution agreement, direct model predictions in both orientations, and edit effects. The natural parent sequence and hash remain on every variant row. DA statistics and observed target/off-target specificity describe the natural parent locus; they are not measured variant properties. Prediction effects describe a 500bp insert within its native 2,114bp context, not an isolated reporter construct. Missing off-target variant predictions remain explicit.

Rebuild from repository root:

```bash
.venv/bin/python scripts/assemble_phase5_library.py --output-dir reports/phase5/library_draft_new
```

The pinned motif-scan config references the checked-in draft and FASTA, so rebuilding into a new directory requires updating the config hashes deliberately. No account or private cluster directory is stored in this config or Slurm script.

The next CPU job scans all sequences with both-strand FIMO against the pinned upstream motif database. It uses nominal site p < 1e-4, the database background, and text output. Text mode does not supply q-values; matches are descriptive rather than FDR-controlled binding calls. See the [FIMO documentation](https://meme-suite.org/meme/doc/fimo.html?man_type=cmd) and [output definition](https://meme-suite.org/meme/doc/fimo-output-format.html).

Outputs retain every motif ID, strand, p-value, position and edited-base overlap. Coordinates convert FIMO's 1-based inclusive positions to 0-based half-open intervals. Parent-versus-variant gained/lost sites refer only to crossing the same nominal threshold; any changed match must overlap the substituted base. The database contains enzyme motifs and redundant TF-family motifs, which require scientific review. A named motif match is not evidence that its named TF binds the candidate or that the edit acts through that motif. Site-level attribution overlap is still pending; existing edit nomination count/profile attribution deltas remain separate evidence.

Submit after pulling the committed code on Great Lakes, with the account and scratch root supplied manually:

```bash
sbatch --account="$SLURM_ACCOUNT" hpc/slurm/scan_phase5_library_motifs.sbatch
```

`PHASE3_SCRATCH_ROOT` must already be exported. The job uses one CPU, 8GB memory, four hours maximum and no GPU. It verifies input hashes, refuses existing output directories, and saves raw FIMO output, annotated hits, the integrated sequence table, runtime/container provenance and output checksums under `phase5_library_motifs` in scratch.

Before synthesis, review motif evidence and flagged natural loci, decide the final parent/variant/control subset, and specify reporter insert length, cloning constraints, adapters and barcodes. This scan does not establish variant off-target specificity or replace experimental validation.
