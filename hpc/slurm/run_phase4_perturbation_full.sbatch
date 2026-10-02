#!/usr/bin/env bash
#SBATCH --job-name=mg_perturb_full
#SBATCH --partition=spgpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus=a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
attribution_source="$scratch_root/phase4_attribution/dual_orientation_attributions.h5"
output_dir="$scratch_root/phase4_perturbation_full"
work_dir="data/intermediate/phase4_perturbation_full_job_${SLURM_JOB_ID}"
model="models/chrombpnet/fold_0/seed_42_neg1/models/muller_chrombpnet_nobias.h5"
candidate_qc="reports/phase4/attribution/attribution_qc.tsv.gz"
parents="reports/phase4/parent_scoring/parent_predictions.tsv.gz"

cd "$project_root"
module load singularity
mkdir -p logs/slurm
for path in "$attribution_source" "$model" "$candidate_qc" "$parents" config/phase4_perturbation_full.json; do
  [[ -s "$path" ]] || { echo "Missing perturbation input: $path" >&2; exit 1; }
done
[[ ! -e "$output_dir" ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 1; }
[[ ! -e "$work_dir" ]] || { echo "Refusing to overwrite: $work_dir" >&2; exit 1; }
[[ "$(sha256sum "$attribution_source" | awk '{print $1}')" == "78d7a6e6c5b875b00de33449131cdd30bbd172e68b791408b4ac114fba4e1b95" ]]
[[ "$(sha256sum "$model" | awk '{print $1}')" == "687befaa442fccae8896326acc884f9578cba5c5d2dd92412aa97f437384e80f" ]]

mkdir -p "$work_dir" "$output_dir.runtime"
cp "$attribution_source" "$work_dir/dual_orientation_attributions.h5"
{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unset}"
  echo "node=$HOSTNAME"
  singularity --version
  nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
  sha256sum containers/chrombpnet-upstream-20260914.sif "$attribution_source" "$model" "$candidate_qc" "$parents"
} > "$output_dir.runtime/runtime.txt"

hpc/run_chrombpnet_container.sh --gpu \
  python scripts/run_phase4_perturbation_pilot.py \
  --candidate-qc "/workspace/$candidate_qc" \
  --parents "/workspace/$parents" \
  --attributions "/workspace/$work_dir/dual_orientation_attributions.h5" \
  --model "/workspace/$model" \
  --config /workspace/config/phase4_perturbation_full.json \
  --output-dir "/workspace/$work_dir/output"

mv "$work_dir/output" "$output_dir"
mv "$output_dir.runtime/runtime.txt" "$output_dir/runtime.txt"
rm "$work_dir/dual_orientation_attributions.h5"
rmdir "$work_dir" "$output_dir.runtime"
(cd "$output_dir" && \
  find . -maxdepth 1 -type f ! -name output_checksums.sha256 -print0 \
    | sort -z | xargs -0 sha256sum > output_checksums.sha256)

echo "Full single-base scan completed. Review coverage and specificity limitations before library selection."
