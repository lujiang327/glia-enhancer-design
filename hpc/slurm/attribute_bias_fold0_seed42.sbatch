#!/usr/bin/env bash
#SBATCH --job-name=mg_biasattr_f0_s42
#SBATCH --partition=spgpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus=a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

# Supply the private account only at submission.
project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
cd "$project_root"
module load singularity

run_name="fold_0/threshold_0.5"
prep_dir="models/bias/$run_name/prep"
train_dir="models/bias/$run_name/seed_42"
numerical_dir="$train_dir/qc_prediction"
attribution_dir="$train_dir/qc_attribution"
model="$train_dir/models/muller_bias.h5"
peaks="$prep_dir/auxiliary/muller_filtered.bias_peaks.bed"

required=(
  "$model"
  "$peaks"
  "$numerical_dir/logs/numerical_qc_summary.json"
)
for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "Missing attribution input: $path" >&2; exit 1; }
done
if ! grep -Eq '"status": "(PASS_NUMERICAL_QC_MOTIF_QC_PENDING|CAUTION_NUMERICAL_QC_MOTIF_QC_REQUIRED)"' \
  "$numerical_dir/logs/numerical_qc_summary.json"; then
  echo "Numerical bias QC did not pass or enter the documented caution range" >&2
  exit 1
fi
if [[ -e "$attribution_dir" ]]; then
  echo "Refusing to overwrite existing output: $attribution_dir" >&2
  exit 1
fi
mkdir -p "$attribution_dir/scores" "$attribution_dir/logs" logs/slurm

hpc/run_chrombpnet_container.sh \
  python scripts/select_bias_motif_regions.py \
  --regions "/workspace/$peaks" \
  --fold /workspace/config/splits/fold_0.json \
  --split test \
  --number 30000 \
  --seed 1234 \
  --output "/workspace/$attribution_dir/test_peaks_30000.bed" \
  --manifest "/workspace/$attribution_dir/logs/region_selection.json"

{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unset}"
  echo "node=$HOSTNAME"
  singularity --version
  nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
  sha256sum containers/chrombpnet-upstream-20260914.sif
  sha256sum "${required[@]}" "$attribution_dir/test_peaks_30000.bed"
} > "$attribution_dir/logs/runtime.txt"

hpc/run_chrombpnet_container.sh --gpu \
  python -m chrombpnet.evaluation.interpret.interpret \
  -g /workspace/data/raw/reference/hg38.canonical.fa \
  -r "/workspace/$attribution_dir/test_peaks_30000.bed" \
  -m "/workspace/$model" \
  -o "/workspace/$attribution_dir/scores/muller_bias" \
  -p counts profile

find "$attribution_dir" -type f ! -name output_checksums.sha256 -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > "$attribution_dir/logs/output_checksums.sha256"

echo "Bias attribution completed. TF-MoDISco and manual motif-leakage review remain pending."
