#!/usr/bin/env bash
#SBATCH --job-name=mg_cbdonor_f0_s42_n1
#SBATCH --partition=spgpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus=a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --array=0-3%2
#SBATCH --output=logs/slurm/%x-%A_%a.out
#SBATCH --error=logs/slurm/%x-%A_%a.err

set -euo pipefail

donors=(LGS1 LGS2 LGS3 LVG1)
donor="${donors[${SLURM_ARRAY_TASK_ID:?Missing array task ID}]}"
project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
cd "$project_root"
module load singularity

prep_dir="models/chrombpnet/fold_0/prep"
train_dir="models/chrombpnet/fold_0/seed_42_neg1"
qc_root="$train_dir/qc_donor_sensitivity"
qc_dir="$qc_root/$donor"
peaks="$prep_dir/auxiliary/muller_filtered.peaks.bed"
nonpeaks="$prep_dir/auxiliary/muller_filtered.nonpeaks.bed"
model="$train_dir/models/muller_chrombpnet.h5"
scaled_bias="$prep_dir/auxiliary/muller_bias_model_scaled.h5"
donor_track="data/intermediate/training/${donor}.insertions.bw"
required=("$peaks" "$nonpeaks" "$model" "$scaled_bias" "$donor_track")
for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "Missing donor-QC input: $path" >&2; exit 1; }
done
if [[ -e "$qc_dir" ]]; then
  echo "Refusing to overwrite existing donor output: $qc_dir" >&2
  exit 1
fi
mkdir -p "$qc_dir/evaluation" "$qc_dir/logs" logs/slurm

{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unset}"
  echo "slurm_array_task_id=${SLURM_ARRAY_TASK_ID:-unset}"
  echo "donor=$donor"
  echo "node=$HOSTNAME"
  singularity --version
  nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
  sha256sum containers/chrombpnet-upstream-20260914.sif
  sha256sum "${required[@]}"
} > "$qc_dir/logs/runtime.txt"

for label in chrombpnet scaled_bias; do
  if [[ "$label" == chrombpnet ]]; then checkpoint="$model"; else checkpoint="$scaled_bias"; fi
  hpc/run_chrombpnet_container.sh --gpu \
    python -m chrombpnet.training.predict \
    -g /workspace/data/raw/reference/hg38.canonical.fa \
    -b "/workspace/$donor_track" \
    -p "/workspace/$peaks" \
    -n "/workspace/$nonpeaks" \
    -o "/workspace/$qc_dir/evaluation/muller_${label}" \
    -fl /workspace/config/splits/fold_0.json \
    -m "/workspace/$checkpoint" \
    -bs 256 -s 42 -il 2114 -ol 1000
done

for path in \
  "$qc_dir/evaluation/muller_chrombpnet_metrics.json" \
  "$qc_dir/evaluation/muller_scaled_bias_metrics.json"; do
  [[ -s "$path" ]] || { echo "Missing donor-QC output: $path" >&2; exit 1; }
done

find "$qc_dir" -type f ! -name output_checksums.sha256 -print0 \
  | sort -z | xargs -0 sha256sum > "$qc_dir/logs/output_checksums.sha256"

echo "Donor $donor held-out prediction QC completed. Aggregate donor review remains required."
