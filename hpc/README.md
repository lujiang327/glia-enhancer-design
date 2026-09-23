# Great Lakes execution

Cluster scripts contain no private account or project path. Submit from the
repository root and provide the account on the command line. Large inputs,
containers, Slurm logs and model outputs stay outside Git.

The checksum-pinned runtime is described in
`config/greatlakes_runtime.lock.json`. `run_chrombpnet_container.sh` verifies the
SIF checksum, activates the pinned ChromBPNet source tree, and applies the
explicit NVIDIA loader path required by the Great Lakes Singularity installation.

## Bias preparation gate

Prepare fold 0 at the initial ATAC bias threshold of 0.5 on a CPU node:

```bash
mkdir -p logs/slurm
sbatch --account=YOUR_ACCOUNT hpc/slurm/prepare_bias_fold0.sbatch
```

This uses the already validated +4/-4 pooled insertion bigWig. It creates a Tn5
PWM diagnostic, filters bias-training regions, and determines the count-loss
weight. It performs no model training. Inspect the Slurm log, generated parameter
TSVs, retained nonpeak count and PWM before submitting a GPU training job.

Output is written under `models/bias/fold_0/threshold_0.5/prep`. The job refuses
to overwrite that directory. Move a failed partial output aside with a descriptive
name before rerunning so its evidence is retained.

## First bias-training candidate

After reviewing the threshold, retained-region count and Tn5 PWM, submit fold 0,
seed 42 on one A40:

```bash
mkdir -p logs/slurm
sbatch --account=YOUR_ACCOUNT hpc/slurm/train_bias_fold0_seed42.sbatch
```

This trains for at most 50 epochs with early-stopping patience 5. It uses the
prepared +4/-4 bigWig directly, avoiding repeated fragment conversion and GPU
time spent sorting reads. It records finite epoch losses and the best validation
epoch. Completion produces a candidate model only; bias prediction, attribution,
motif leakage checks and peak-transfer QC remain required before ChromBPNet
training.

## Held-out numerical bias QC

After successful candidate training, evaluate peaks and nonpeaks on the fold 0
test chromosomes (`chr1`, `chr3`, and `chr6`):

```bash
mkdir -p logs/slurm
sbatch --account=YOUR_ACCOUNT hpc/slurm/evaluate_bias_fold0_seed42.sbatch
```

This writes separate held-out count correlations and profile JSD metrics for
peaks and nonpeaks. A positive nonpeak count correlation and peak correlation
above -0.3 pass the preliminary screen. Peak correlation from -0.5 through
-0.3 is retained as a caution; at or below -0.5 fails. Motif-leakage QC is a
separate required gate.

## Bias attribution and motif-leakage gate

After a numerical pass, compute counts and profile contribution scores on a
deterministic sample of 30,000 held-out peaks:

```bash
sbatch --account=YOUR_ACCOUNT hpc/slurm/attribute_bias_fold0_seed42.sbatch
```

After attribution completes, run TF-MoDISco on a CPU node, preferably with an
`afterok` dependency on the attribution job:

```bash
sbatch --account=YOUR_ACCOUNT hpc/slurm/modisco_bias_fold0_seed42.sbatch
```

The two HTML reports must be reviewed before accepting the bias model. Dominant
Tn5/enzyme-bias patterns or unstructured repeats are expected. An obvious
transcription-factor motif in the leading patterns blocks use of this bias model
until the threshold/model choice is revisited.

## Full ChromBPNet preparation gate

After the bias model passes numerical and motif review, prepare fold 0 full-model
regions, count-loss weight, and a depth-scaled frozen bias model:

```bash
sbatch --account=YOUR_ACCOUNT hpc/slurm/prepare_chrombpnet_fold0.sbatch
```

This job uses the GPU only to scale the accepted bias model; it does not train
the full model. It uses a fixed seed for the upstream nonpeak subsample, retains
500-bp training jitter, and writes all preparation artifacts under
`models/chrombpnet/fold_0/prep`. Review the retained region counts, thresholds,
scaled-model checksum, and parameters before submitting full-model training.

## First full ChromBPNet candidate

After reviewing the fold 0 preparation summary, train seed 42 on one A40:

```bash
sbatch --account=YOUR_ACCOUNT hpc/slurm/train_chrombpnet_fold0_seed42.sbatch
```

The job uses genome-wide prepared peaks plus matched backgrounds, the frozen
depth-scaled bias model, 500-bp jitter, and chromosome-disjoint train/validation
sets. It produces both combined and no-bias checkpoints. Completion is not a
model pass: held-out profile/count performance, additional seeds, donor
sensitivity, and attribution/motif recovery remain required.

## Initial held-out full-model numerical comparison

After successful seed-42 training, evaluate the combined model and its scaled
bias baseline on the identical fold-0 test regions (`chr1`, `chr3`, and `chr6`):

```bash
sbatch --account=YOUR_ACCOUNT hpc/slurm/evaluate_chrombpnet_fold0_seed42.sbatch
```

The report records peak, nonpeak, and combined count correlations/MSE and profile
Jensen-Shannon distances. It expresses every model-minus-baseline comparison so
that positive values favor the full model. Peak metric regressions and nonpeak
count-calibration regressions are retained as explicit cautions. This is an
initial aggregate screen;
it does not pass the model gate without per-chromosome, depth, donor, seed,
attribution, and motif assessment.

## Background-calibration control

If the initial candidate improves held-out peaks but broadly inflates nonpeak
counts, train a controlled candidate that changes only
`negative_sampling_ratio` from 0.1 to 1.0:

```bash
sbatch --account=YOUR_ACCOUNT hpc/slurm/train_chrombpnet_fold0_seed42_neg1.sbatch
```

The job writes to `models/chrombpnet/fold_0/seed_42_neg1` and preserves the
original checkpoint. After it converges, evaluate it on the same test regions:

```bash
sbatch --account=YOUR_ACCOUNT hpc/slurm/evaluate_chrombpnet_fold0_seed42_neg1.sbatch
```

Select between the two candidates using their matched held-out peak and nonpeak
metrics. Do not choose on training loss alone.

## Provisional-model attribution and motif review

For the ratio-1 calibration candidate, compute contribution scores from the
no-bias biological checkpoint on a deterministic sample of 30,000 held-out
peaks:

```bash
sbatch --account=YOUR_ACCOUNT hpc/slurm/attribute_chrombpnet_fold0_seed42_neg1.sbatch
```

Then run TF-MoDISco on a CPU node, preferably using an `afterok` dependency:

```bash
sbatch --account=YOUR_ACCOUNT hpc/slurm/modisco_chrombpnet_fold0_seed42_neg1.sbatch
```

Review profile and count reports for coherent Müller-glia regulatory motifs,
seqlet support, diffuse GC/repeat patterns, and residual Tn5 motifs. Motif names
identify families rather than exact expressed TFs; paired RNA provides supporting
evidence. Attribution review is part of model validation and does not open the
differential-accessibility gate by itself.

## Balanced-model seed reproducibility

After the seed-42 balanced candidate passes attribution and motif review, train
two independent replicas and submit their held-out evaluations as a dependent
array:

```bash
TRAIN_JOB=$(sbatch --parsable --account=YOUR_ACCOUNT \
  hpc/slurm/train_chrombpnet_fold0_neg1_seeds.sbatch)
TRAIN_JOB=${TRAIN_JOB%%;*}
sbatch --account=YOUR_ACCOUNT --dependency="afterok:${TRAIN_JOB}" \
  hpc/slurm/evaluate_chrombpnet_fold0_neg1_seeds.sbatch
```

Array tasks 0 and 1 use seeds 123 and 456. They write separate checkpoints and
QC under `seed_123_neg1` and `seed_456_neg1`. Compare both with seed 42 before
starting donor-sensitivity work. Differential accessibility remains blocked
until both reproducibility gates are reviewed.

## Donor-sensitivity gate

After seed reproducibility passes, evaluate the pooled seed-42 balanced model
against each donor's own insertion track on the same held-out chromosomes and
regions. Submit the four-donor GPU array, then submit the CPU summarizer with an
`afterok` dependency:

```bash
DONOR_JOB=$(sbatch --parsable --account=YOUR_ACCOUNT \
  hpc/slurm/evaluate_chrombpnet_fold0_seed42_neg1_donors.sbatch)
DONOR_JOB=${DONOR_JOB%%;*}
SUMMARY_JOB=$(sbatch --parsable --account=YOUR_ACCOUNT \
  --dependency="afterok:${DONOR_JOB}" \
  hpc/slurm/summarize_chrombpnet_fold0_seed42_neg1_donors.sbatch)
SUMMARY_JOB=${SUMMARY_JOB%%;*}
printf 'DONOR_JOB=%s\nSUMMARY_JOB=%s\n' "$DONOR_JOB" "$SUMMARY_JOB"
```

Array tasks 0 through 3 correspond to LGS1, LGS2, LGS3, and LVG1. Two tasks
may run concurrently. The evaluation compares the combined model and matched
scaled-bias baseline with each donor's observed insertions. Count correlations
and profile metrics are comparable across donors. Absolute count MSE is saved
but excluded from the donor decision because the checkpoint's count scale was
learned from the pooled library. This is donor-resolved evaluation of a pooled
model, not leave-one-donor-out retraining.
