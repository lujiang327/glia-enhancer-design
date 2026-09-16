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
