# Model-first glial regulatory-element workflow

The current deliverable is a data-feasibility assessment, not a trained model.
QC and differential accessibility are conceptually separate analyses. Project
execution nevertheless places training and validation before all differential
accessibility, candidate selection, optimization, and experimental-library work.

## Phase 1: auditable QC

Use the hg38 GSE196235 published Müller glia fragment pool first. Its physical GEO
attachment is GSM5866073 despite the series description referring to rep1. Never
interpret the attachment sample as the donor of every fragment. The published
RNA counts retain sample prefixes; the fragment file does not. Match barcodes
only when exactly one RNA sample contains the bare barcode. Exclude every
ambiguous barcode from the donor-resolved working cohort, without assigning a
new cell type. Keep the original complete fragment file and the audit of all
excluded barcodes. Both eyes from the same person are one biological donor.

The sample-prefix mapping comes from the archived authors' scRNA code; GEO maps
the sample names to accessions. RNA label evidence is membership in the authors'
Müller glia fragment pool plus an unambiguous paired RNA barcode. Check GLUL and
RLBP1 with SLC1A3/SOX9 and other retinal markers. Mixed neuronal RNA can reflect
ambient RNA or residual doublets; expression alone does not justify relabeling.

Count each fragment row once, irrespective of read-support column 5. Do not
deduplicate across biological cells at identical genomic coordinates. Exclude
chrM, chrY and alternative contigs in the current chr1–22,X policy; remove any
fragment overlapping hg38 blacklist v2, and out-of-bounds coordinates. These are
retained processed records, not a fresh verification of alignment MAPQ or library
complexity, which require original alignments. Keep both full-pool and cohort
counts explicitly distinguished. Retain donor pseudobulks for future use.

TSS scores computed here are pooled RefGene insertion-enrichment proxies, not
ArchR per-cell scores. FRiP uses one count per fragment with any overlap against
the diagnostic peak set and the same genomic-filtered denominator. Diagnostic
peak calls use MACS3 BEDPE, hs effective genome, q=0.01, keep-dup all; no DA test
is involved. Such peaks are not yet the final relaxed insertion-based training
peak universe. Compare donors using blacklist-free 50-kb midpoint counts and
explicitly label these coarse QC bins, not accessible peaks. Four donors cannot
cleanly disentangle donor biology from technical batch.

## Phase 2: ChromBPNet preparation and validation

Use a Linux NVIDIA/CUDA environment with the pinned ChromBPNet commit in
config/training_plan.json. Confirm actual GPU detection, available RAM/scratch,
TensorFlow/CUDA compatibility, and container digest before submitting jobs.
Do not silently install an incompatible Apple-Silicon training stack. Save
environment exports, package versions, commands, stdout/stderr and job IDs.

Obtain and checksum an hg38 analysis-set FASTA matching the coordinates and
chromosome sizes. Use the official fragment-to-insertion preprocessing and save
its automatic enzyme-shift estimation. The source 10x fragments already contain
Tn5 offsets: do not blindly add +4/-5 again. The current TSS proxy uses start and
end-1; the final model track must follow the pinned implementation's documented
convention and shift diagnostics. Keep molecule counts unweighted by read support.

Recall genome-wide relaxed accessible regions for the working cohort using the
official ChromBPNet preprocessing recipe and retain summit-centered narrowPeak
records. Exclude blacklist, invalid/reference-N and boundary windows, including
the full input width plus jitter. Generate GC-matched genomic nonpeak backgrounds
excluding accessible regions and blacklist; freeze all regions and their hashes.
Never substitute the later ~1,000 candidate regions for this training universe.

Chromosome fold JSONs in config/splits reuse the published retina partitions,
restricted to chr1–22,X (see make_splits.py). The same fold applies
to peaks, negatives, bias fitting, tuning and evaluation. No interval, reverse
complement or jittered version can cross folds; assess segmental-duplication
homology separately because chromosome splits alone do not remove all homologous
sequence. Start with fold 0 across seeds 42/123/456. Use only training and validation
chromosomes to choose hyperparameters; use test chromosomes once decisions are
frozen. Rotate all five folds for the final model ensemble. Leave one biological
donor out in additional fits and compare against that donor on held-out chromosomes.

Fit a multiome-compatible Tn5 bias model on nonpeaks (initial threshold factor
0.5). Review GC distributions, count/profile QC, and motifs. A bias model learning
cell-type TF motifs is a failure; a TF model retaining Tn5 bias is also a failure.
Tune bias threshold on training/validation data per official guidance; archive
each attempt instead of overwriting failed runs.

Prespecify the evaluation table before training: epoch losses, early-stop epoch,
finite gradients/losses, peak-only and peak-plus-background count Pearson and
Spearman on log1p insertion counts, count calibration, multinomial profile loss,
Jensen–Shannon profile distance with a fixed definition/smoothing, and results by
chromosome, depth stratum and donor. Compare profiles to bias-only and flat-profile
baselines; compare counts to training-derived simple baselines. Bootstrap genomic
blocks for uncertainty and compare to donor/split-half noise ceilings. Do not
declare success from an arbitrary universal depth or correlation cutoff.

Require reproducible improvement over relevant baselines, no collapse/large seed
dependence, sensible held-out donor behavior, and biologically interpretable
attribution/motifs. Run TF-MoDISco and motif-family matching; use paired RNA to
support which family members are expressed. Do not claim an exact TF identity
from a shared motif. Inspect positive and negative attribution and whether it is
stable across seeds. If any check fails, save the failed result and stop. Only a
documented validated checkpoint with completed evidence can pass model_gate.py.
The gate records scientific review; it does not replace review with file existence.

## Phases 3–5: disabled until a usable model exists

After validation, perform donor-replicated pseudobulk DA with a design accounting
for paired donor comparisons and relevant covariates. Use Müller glia vs rest and
important pairwise off-target contrasts (including astrocytes, rods, cones, bipolar,
amacrine, ganglion, horizontal and microglial classes where adequate). Handle
missing donor/cell-type strata explicitly; cells and eyes are not independent
biological replicates. Define how the rest mixture is weighted. Rank ~1,000 open
regions by effect size, FDR, target/off-target accessibility and reproducibility.
Annotate promoter, intronic, distal and other genomic classes without requiring
enhancer confirmation. DA statistics must be computed from the appropriate
donor count matrix, not ChromBPNet predictions.

Attribute and perturb candidates using the validated model. Keep the complete
reference context (2114 bp by default) distinct from the eventual assay insert;
define insert length before library manufacture. Preserve every natural parent,
genome/reference hash, strand, edit list and parent identifier. Compare sequence
variants across seeds/folds and retain uncertainty. Maximizing Müller glia signal
alone is not a specificity claim. Use calibrated off-target model predictions or
measured accessibility where available; otherwise leave specificity unavailable.

The future library schema includes: candidate_id, parent_id, genome_build, chrom,
start, end, strand, genomic_class, parent_sequence, variant_sequence, edits,
DA_contrast, log2FC, FDR, target_accessibility, off_target_accessibility,
donor_reproducibility, model_id, model_prediction, attribution/motif_summary,
predicted_delta_vs_parent, specificity_score, specificity_method, uncertainty.
No fabricated sequence or placeholder statistics are inserted into a deliverable.

## Reuse for human brain astrocytes

Only after retinal validation, add Allen/CELLxGENE dataset adapters with explicit
fragment-to-annotation barcode mapping, adult tissue/region/donor metadata and
genome/assay checks. RNA-only CELLxGENE assets cannot substitute for ATAC fragments.
Reuse the gate, QC, donor construction, folds and validation; replace retinal
marker/off-target panels with astrocyte/brain-appropriate ones. Do not transfer
retinal donor identities, performance conclusions or training checkpoints as
evidence that a brain model works.
