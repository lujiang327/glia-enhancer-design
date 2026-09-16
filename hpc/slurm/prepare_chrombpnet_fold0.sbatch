#!/usr/bin/env bash
#SBATCH --job-name=mg_cbprep_f0
#SBATCH --partition=spgpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus=a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

# Supply the private account only at submission.
project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
cd "$project_root"
module load singularity

prep_dir="models/chrombpnet/fold_0/prep"
bias_model="models/bias/fold_0/threshold_0.5/seed_42/models/muller_bias.h5"
required=(
  "$bias_model"
  data/raw/reference/hg38.canonical.fa
  data/raw/reference/hg38.canonical.fa.fai
  data/intermediate/training/muller.insertions.bw
  data/intermediate/training/muller.training_peaks.narrowPeak
  data/intermediate/training/fold_0.negatives.bed
  config/splits/fold_0.json
)
for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "Missing ChromBPNet-preparation input: $path" >&2; exit 1; }
done
if [[ -e "$prep_dir" ]]; then
  echo "Refusing to overwrite existing output: $prep_dir" >&2
  exit 1
fi
mkdir -p "$prep_dir/auxiliary" "$prep_dir/logs" logs/slurm

{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unset}"
  echo "node=$HOSTNAME"
  singularity --version
  nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
  sha256sum containers/chrombpnet-upstream-20260914.sif
  sha256sum "${required[@]}"
} > "$prep_dir/logs/runtime.txt"

# This performs no full-model training. It filters train/validation outliers,
# fixes the count-loss weight, and scales the accepted bias model to the pooled
# Müller-glia depth. The wrapper supplies the RNG seed missing from upstream.
hpc/run_chrombpnet_container.sh --gpu \
  python scripts/run_chrombpnet_hyperparams.py \
  -g /workspace/data/raw/reference/hg38.canonical.fa \
  -i /workspace/data/intermediate/training/muller.insertions.bw \
  -p /workspace/data/intermediate/training/muller.training_peaks.narrowPeak \
  -n /workspace/data/intermediate/training/fold_0.negatives.bed \
  -sr 0.1 \
  -oth 0.9999 \
  -j 500 \
  -fl /workspace/config/splits/fold_0.json \
  -il 2114 \
  -ol 1000 \
  -fil 512 \
  -dil 8 \
  -b "/workspace/$bias_model" \
  -op "/workspace/$prep_dir/auxiliary/muller_" \
  --sampling-seed 42

outputs=(
  "$prep_dir/auxiliary/muller_filtered.peaks.bed"
  "$prep_dir/auxiliary/muller_filtered.nonpeaks.bed"
  "$prep_dir/auxiliary/muller_bias_model_scaled.h5"
  "$prep_dir/auxiliary/muller_chrombpnet_data_params.tsv"
  "$prep_dir/auxiliary/muller_chrombpnet_model_params.tsv"
)
for path in "${outputs[@]}"; do
  [[ -s "$path" ]] || { echo "Missing ChromBPNet-preparation output: $path" >&2; exit 1; }
done
{
  echo "filtered_peaks=$(wc -l < "$prep_dir/auxiliary/muller_filtered.peaks.bed")"
  echo "filtered_nonpeaks=$(wc -l < "$prep_dir/auxiliary/muller_filtered.nonpeaks.bed")"
  cat "$prep_dir/auxiliary/muller_chrombpnet_data_params.tsv"
  cat "$prep_dir/auxiliary/muller_chrombpnet_model_params.tsv"
} > "$prep_dir/logs/preparation_summary.txt"
sha256sum "${outputs[@]}" > "$prep_dir/logs/output_checksums.sha256"

echo "Full-model preparation completed. Inspect parameters before training ChromBPNet."
