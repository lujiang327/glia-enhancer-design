#!/usr/bin/env bash
#SBATCH --job-name=mg_biasqc_f0_s42
#SBATCH --partition=spgpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus=a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

# Supply the private account only at submission:
# sbatch --account=YOUR_ACCOUNT hpc/slurm/evaluate_bias_fold0_seed42.sbatch

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
cd "$project_root"
module load singularity

run_name="fold_0/threshold_0.5"
prep_dir="models/bias/$run_name/prep"
train_dir="models/bias/$run_name/seed_42"
qc_dir="$train_dir/qc_prediction"

required=(
  "$train_dir/models/muller_bias.h5"
  "$prep_dir/auxiliary/muller_filtered.bias_nonpeaks.bed"
  "$prep_dir/auxiliary/muller_filtered.bias_peaks.bed"
)
for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "Missing bias-QC input: $path" >&2; exit 1; }
done
if [[ -e "$qc_dir" ]]; then
  echo "Refusing to overwrite existing output: $qc_dir" >&2
  exit 1
fi
mkdir -p "$qc_dir/evaluation" "$qc_dir/logs" logs/slurm

{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unset}"
  echo "node=$HOSTNAME"
  singularity --version
  nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
  sha256sum containers/chrombpnet-upstream-20260914.sif
  sha256sum "${required[@]}"
} > "$qc_dir/logs/runtime.txt"

# ChromBPNet's predictor restricts evaluation to the test chromosomes in the
# fold definition. Supplying both region sets yields separate peak/nonpeak
# metrics as well as combined metrics; no training or validation regions enter.
hpc/run_chrombpnet_container.sh --gpu \
  python -m chrombpnet.training.predict \
  -g /workspace/data/raw/reference/hg38.canonical.fa \
  -b /workspace/data/intermediate/training/muller.insertions.bw \
  -p /workspace/models/bias/$run_name/prep/auxiliary/muller_filtered.bias_peaks.bed \
  -n /workspace/models/bias/$run_name/prep/auxiliary/muller_filtered.bias_nonpeaks.bed \
  -o /workspace/models/bias/$run_name/seed_42/qc_prediction/evaluation/muller_bias \
  -fl /workspace/config/splits/fold_0.json \
  -m /workspace/models/bias/$run_name/seed_42/models/muller_bias.h5 \
  -bs 256 \
  -s 42 \
  -il 2114 \
  -ol 1000

hpc/run_chrombpnet_container.sh \
  python scripts/summarize_bias_prediction_qc.py \
  --metrics /workspace/models/bias/$run_name/seed_42/qc_prediction/evaluation/muller_bias_metrics.json \
  --output /workspace/models/bias/$run_name/seed_42/qc_prediction/logs/numerical_qc_summary.json

find "$qc_dir" -type f ! -name output_checksums.sha256 -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > "$qc_dir/logs/output_checksums.sha256"

echo "Held-out numerical bias QC completed. Attribution/motif QC remains pending."
