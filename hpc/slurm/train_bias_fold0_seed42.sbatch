#!/usr/bin/env bash
#SBATCH --job-name=mg_bias_f0_s42_b05
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

# Supply the private account only at submission:
# sbatch --account=YOUR_ACCOUNT hpc/slurm/train_bias_fold0_seed42.sbatch

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
cd "$project_root"
module load singularity

run_name="fold_0/threshold_0.5"
prep_dir="models/bias/$run_name/prep"
train_dir="models/bias/$run_name/seed_42"

required=(
  "$prep_dir/logs/muller_bias_data_params.tsv"
  "$prep_dir/logs/muller_bias_model_params.tsv"
  "$prep_dir/auxiliary/muller_filtered.bias_nonpeaks.bed"
  "$prep_dir/auxiliary/muller_filtered.bias_peaks.bed"
)
for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "Missing bias-preparation output: $path" >&2; exit 1; }
done
if [[ -e "$train_dir" ]]; then
  echo "Refusing to overwrite existing output: $train_dir" >&2
  exit 1
fi
mkdir -p "$train_dir/models" "$train_dir/logs" logs/slurm

{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unset}"
  echo "node=$HOSTNAME"
  singularity --version
  nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
  sha256sum containers/chrombpnet-upstream-20260914.sif
  sha256sum "${required[@]}"
} > "$train_dir/logs/runtime.txt"

# This is the first full bias-training candidate. The upstream defaults are
# explicit here: 50 maximum epochs, early stopping patience 5, batch size 64,
# learning rate 0.001, 128 filters and four dilation layers (recorded in the
# preparation parameter file). QC and motif acceptance run only after this job.
hpc/run_chrombpnet_container.sh --gpu \
  python -m chrombpnet.training.train \
  -g /workspace/data/raw/reference/hg38.canonical.fa \
  -b /workspace/data/intermediate/training/muller.insertions.bw \
  -p None \
  -n /workspace/models/bias/$run_name/prep/auxiliary/muller_filtered.bias_nonpeaks.bed \
  -o /workspace/models/bias/$run_name/seed_42/models/muller_bias \
  -fl /workspace/config/splits/fold_0.json \
  -e 50 \
  -es 5 \
  -bs 64 \
  -l 0.001 \
  -pf /workspace/models/bias/$run_name/prep/logs/muller_bias_model_params.tsv \
  -s 42 \
  -a /workspace/data/intermediate/vendor/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c/chrombpnet/training/models/bpnet_model.py \
  --trackables \
    logcount_predictions_loss \
    loss \
    logits_profile_predictions_loss \
    val_logcount_predictions_loss \
    val_loss \
    val_logits_profile_predictions_loss

mv "$train_dir/models/muller_bias.log" "$train_dir/logs/"
mv "$train_dir/models/muller_bias.log.batch" "$train_dir/logs/"
mv "$train_dir/models/muller_bias.args.json" "$train_dir/logs/"

hpc/run_chrombpnet_container.sh --gpu \
  python scripts/summarize_keras_history.py \
  --history /workspace/models/bias/$run_name/seed_42/logs/muller_bias.log \
  --output /workspace/models/bias/$run_name/seed_42/logs/training_summary.json

sha256sum \
  "$train_dir/models/muller_bias.h5" \
  "$train_dir/logs/muller_bias.log" \
  "$train_dir/logs/muller_bias.log.batch" \
  "$train_dir/logs/muller_bias.args.json" \
  "$train_dir/logs/training_summary.json" \
  > "$train_dir/logs/output_checksums.sha256"

echo "Bias training completed. Model remains a candidate pending bias QC."
