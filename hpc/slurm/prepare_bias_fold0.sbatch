#!/usr/bin/env bash
#SBATCH --job-name=mg_biasprep_f0_b05
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

# Supply the private account only at submission:
# sbatch --account=YOUR_ACCOUNT hpc/slurm/prepare_bias_fold0.sbatch

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
cd "$project_root"
module load singularity

run_name="fold_0/threshold_0.5"
prep_dir="models/bias/$run_name/prep"
if [[ -e "$prep_dir" ]]; then
  echo "Refusing to overwrite existing output: $prep_dir" >&2
  exit 1
fi
mkdir -p "$prep_dir/auxiliary" "$prep_dir/evaluation" "$prep_dir/logs" logs/slurm

{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unset}"
  echo "node=$HOSTNAME"
  singularity --version
  sha256sum containers/chrombpnet-upstream-20260914.sif
} > "$prep_dir/logs/runtime.txt"

# The prepared pooled bigWig is already +4/-4 normalized and independently
# validated against the four donor tracks. Build the shift-PWM diagnostic from
# fold 0's first validation chromosome, then derive the thresholded bias regions
# and model parameters without regenerating the signal from fragments.
hpc/run_chrombpnet_container.sh \
  python -m chrombpnet.helpers.preprocessing.analysis.build_pwm_from_bigwig \
  -i /workspace/data/intermediate/training/muller.insertions.bw \
  -g /workspace/data/raw/reference/hg38.canonical.fa \
  -op /workspace/models/bias/$run_name/prep/evaluation/muller_bw_shift_qc \
  -cr chr8 \
  -c /workspace/data/raw/reference/hg38.canonical.chrom.sizes \
  -pw 24

hpc/run_chrombpnet_container.sh \
  python -m chrombpnet.helpers.hyperparameters.find_bias_hyperparams \
  -g /workspace/data/raw/reference/hg38.canonical.fa \
  -i /workspace/data/intermediate/training/muller.insertions.bw \
  -p /workspace/data/intermediate/training/muller.training_peaks.narrowPeak \
  -n /workspace/data/intermediate/training/fold_0.negatives.bed \
  -b 0.5 \
  -oth 0.9999 \
  -j 0 \
  -fl /workspace/config/splits/fold_0.json \
  -il 2114 \
  -ol 1000 \
  -fil 128 \
  -dil 4 \
  -op /workspace/models/bias/$run_name/prep/auxiliary/muller_

mv "$prep_dir/auxiliary/muller_bias_data_params.tsv" "$prep_dir/logs/"
mv "$prep_dir/auxiliary/muller_bias_model_params.tsv" "$prep_dir/logs/"
sha256sum \
  "$prep_dir/logs/muller_bias_data_params.tsv" \
  "$prep_dir/logs/muller_bias_model_params.tsv" \
  "$prep_dir/auxiliary/muller_filtered.bias_nonpeaks.bed" \
  "$prep_dir/auxiliary/muller_filtered.bias_peaks.bed" \
  > "$prep_dir/logs/output_checksums.sha256"

echo "Bias preparation completed. Inspect thresholds and region counts before GPU training."
