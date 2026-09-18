#!/usr/bin/env bash
#SBATCH --job-name=mg_cbqc_f0_s42
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

# Supply the private account only at submission.
project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
cd "$project_root"
module load singularity

prep_dir="models/chrombpnet/fold_0/prep"
train_dir="models/chrombpnet/fold_0/seed_42"
qc_dir="$train_dir/qc_prediction"
peaks="$prep_dir/auxiliary/muller_filtered.peaks.bed"
nonpeaks="$prep_dir/auxiliary/muller_filtered.nonpeaks.bed"
model="$train_dir/models/muller_chrombpnet.h5"
scaled_bias="$prep_dir/auxiliary/muller_bias_model_scaled.h5"
required=("$peaks" "$nonpeaks" "$model" "$scaled_bias")
for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "Missing full-model QC input: $path" >&2; exit 1; }
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

# Evaluate the combined model and the depth-scaled frozen bias model on exactly
# the same fold-0 test regions. The predictor reads only chr1, chr3, and chr6
# from the supplied peak and nonpeak files according to fold_0.json.
for label in chrombpnet scaled_bias; do
  if [[ "$label" == chrombpnet ]]; then
    checkpoint="$model"
  else
    checkpoint="$scaled_bias"
  fi
  hpc/run_chrombpnet_container.sh --gpu \
    python -m chrombpnet.training.predict \
    -g /workspace/data/raw/reference/hg38.canonical.fa \
    -b /workspace/data/intermediate/training/muller.insertions.bw \
    -p /workspace/$peaks \
    -n /workspace/$nonpeaks \
    -o /workspace/$qc_dir/evaluation/muller_${label} \
    -fl /workspace/config/splits/fold_0.json \
    -m /workspace/$checkpoint \
    -bs 256 \
    -s 42 \
    -il 2114 \
    -ol 1000
done

hpc/run_chrombpnet_container.sh \
  python scripts/compare_chrombpnet_prediction_qc.py \
  --model-metrics /workspace/$qc_dir/evaluation/muller_chrombpnet_metrics.json \
  --bias-metrics /workspace/$qc_dir/evaluation/muller_scaled_bias_metrics.json \
  --output /workspace/$qc_dir/logs/numerical_qc_summary.json

find "$qc_dir" -type f ! -name output_checksums.sha256 -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > "$qc_dir/logs/output_checksums.sha256"

echo "Held-out full-model numerical comparison completed. Scientific review remains required."
