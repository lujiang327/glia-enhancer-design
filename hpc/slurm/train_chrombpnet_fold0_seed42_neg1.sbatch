#!/usr/bin/env bash
#SBATCH --job-name=mg_cb_f0_s42_n1
#SBATCH --partition=spgpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus=a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

# Controlled calibration candidate: identical to seed_42 except that every
# training/validation peak is paired with one nonpeak instead of 0.1 nonpeak.
project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
cd "$project_root"
module load singularity

prep_dir="models/chrombpnet/fold_0/prep"
train_dir="models/chrombpnet/fold_0/seed_42_neg1"
source_params="$prep_dir/auxiliary/muller_chrombpnet_model_params.tsv"
run_params="$train_dir/config/muller_chrombpnet_model_params.tsv"
architecture="data/intermediate/vendor/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c/chrombpnet/training/models/chrombpnet_with_bias_model.py"
required=(
  "$prep_dir/auxiliary/muller_filtered.peaks.bed"
  "$prep_dir/auxiliary/muller_filtered.nonpeaks.bed"
  "$prep_dir/auxiliary/muller_bias_model_scaled.h5"
  "$source_params"
  "$architecture"
)
for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "Missing calibration-training input: $path" >&2; exit 1; }
done
if [[ -e "$train_dir" ]]; then
  echo "Refusing to overwrite existing output: $train_dir" >&2
  exit 1
fi
mkdir -p "$train_dir/config" "$train_dir/models" "$train_dir/logs" logs/slurm

awk 'BEGIN { FS=OFS="\t" }
     $1 == "negative_sampling_ratio" { $2="1.0"; found++ }
     { print }
     END { if (found != 1) exit 2 }' \
  "$source_params" > "$run_params"
grep -Fx $'negative_sampling_ratio\t1.0' "$run_params" >/dev/null

{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unset}"
  echo "node=$HOSTNAME"
  singularity --version
  nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
  sha256sum containers/chrombpnet-upstream-20260914.sif
  sha256sum "${required[@]}" "$run_params"
} > "$train_dir/logs/runtime.txt"

hpc/run_chrombpnet_container.sh --gpu \
  python -m chrombpnet.training.train \
  -g /workspace/data/raw/reference/hg38.canonical.fa \
  -b /workspace/data/intermediate/training/muller.insertions.bw \
  -p /workspace/models/chrombpnet/fold_0/prep/auxiliary/muller_filtered.peaks.bed \
  -n /workspace/models/chrombpnet/fold_0/prep/auxiliary/muller_filtered.nonpeaks.bed \
  -o /workspace/models/chrombpnet/fold_0/seed_42_neg1/models/muller_chrombpnet \
  -fl /workspace/config/splits/fold_0.json \
  -e 50 \
  -es 5 \
  -bs 64 \
  -l 0.001 \
  -pf /workspace/models/chrombpnet/fold_0/seed_42_neg1/config/muller_chrombpnet_model_params.tsv \
  -s 42 \
  -a "/workspace/$architecture" \
  --trackables \
    logcount_predictions_loss \
    loss \
    logits_profile_predictions_loss \
    val_logcount_predictions_loss \
    val_loss \
    val_logits_profile_predictions_loss

outputs=(
  "$train_dir/models/muller_chrombpnet.h5"
  "$train_dir/models/muller_chrombpnet_nobias.h5"
  "$train_dir/models/muller_chrombpnet.log"
  "$train_dir/models/muller_chrombpnet.log.batch"
  "$train_dir/models/muller_chrombpnet.args.json"
)
for path in "${outputs[@]}"; do
  [[ -s "$path" ]] || { echo "Missing calibration-training output: $path" >&2; exit 1; }
done

mv "$train_dir/models/muller_chrombpnet.log" "$train_dir/logs/"
mv "$train_dir/models/muller_chrombpnet.log.batch" "$train_dir/logs/"
mv "$train_dir/models/muller_chrombpnet.args.json" "$train_dir/logs/"

hpc/run_chrombpnet_container.sh \
  python scripts/summarize_keras_history.py \
  --history /workspace/models/chrombpnet/fold_0/seed_42_neg1/logs/muller_chrombpnet.log \
  --output /workspace/models/chrombpnet/fold_0/seed_42_neg1/logs/training_summary.json

find "$train_dir" -type f ! -name output_checksums.sha256 -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > "$train_dir/logs/output_checksums.sha256"

echo "Negative-ratio-1 calibration candidate completed. Held-out comparison remains required."
