#!/usr/bin/env bash
#SBATCH --job-name=mg_parent_attr
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

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
output_dir="$scratch_root/phase4_attribution"
staging_dir="data/intermediate/phase4_attribution_job_${SLURM_JOB_ID}"
candidates="reports/phase4/parent_scoring/parent_predictions.tsv.gz"
genome="data/raw/reference/hg38.canonical.fa"
model="models/chrombpnet/fold_0/seed_42_neg1/models/muller_chrombpnet_nobias.h5"

cd "$project_root"
module load singularity
mkdir -p logs/slurm
for path in "$candidates" "$genome" "$genome.fai" "$model" config/phase4_attribution.json; do
  [[ -s "$path" ]] || { echo "Missing attribution input: $path" >&2; exit 1; }
done
[[ ! -e "$output_dir" ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 1; }
[[ ! -e "$staging_dir" ]] || { echo "Refusing to overwrite staging directory: $staging_dir" >&2; exit 1; }
[[ "$(sha256sum "$genome" | awk '{print $1}')" == "fb96c995cd3a4b40b73873ff33ff885dd32da7e1b5490481aa4cc3a264965cf1" ]]
[[ "$(sha256sum "$model" | awk '{print $1}')" == "687befaa442fccae8896326acc884f9578cba5c5d2dd92412aa97f437384e80f" ]]

mkdir -p "$output_dir.runtime"
{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unset}"
  echo "node=$HOSTNAME"
  singularity --version
  nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
  sha256sum containers/chrombpnet-upstream-20260914.sif "$candidates" "$genome" "$model"
} > "$output_dir.runtime/runtime.txt"

hpc/run_chrombpnet_container.sh --gpu \
  python scripts/attribute_phase4_parent_sequences.py \
  --candidates "/workspace/$candidates" \
  --genome "/workspace/$genome" \
  --model "/workspace/$model" \
  --config /workspace/config/phase4_attribution.json \
  --output-dir "/workspace/$staging_dir"

mv "$staging_dir" "$output_dir"
mv "$output_dir.runtime/runtime.txt" "$output_dir/runtime.txt"
rmdir "$output_dir.runtime"
(cd "$output_dir" && \
  find . -maxdepth 1 -type f ! -name output_checksums.sha256 -print0 \
    | sort -z | xargs -0 sha256sum > output_checksums.sha256)

echo "Dual-orientation attribution completed. Variant generation remains blocked pending attribution QC."
