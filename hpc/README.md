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
