#!/usr/bin/env bash
#SBATCH --job-name=mg_motif_attr
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Set scratch root}"
source_h5="$scratch_root/phase4_attribution/dual_orientation_attributions.h5"
output_dir="$scratch_root/phase5_motif_attribution"
work_dir="data/intermediate/phase5_motif_attribution_job_${SLURM_JOB_ID}"
[[ -s "$source_h5" && ! -e "$output_dir" && ! -e "$work_dir" ]] || { echo 'Missing input or output already exists' >&2; exit 1; }
module load singularity
mkdir -p "$work_dir"
cp "$source_h5" "$work_dir/dual_orientation_attributions.h5"
hpc/run_chrombpnet_container.sh python -u scripts/annotate_phase5_motif_attribution.py \
  --config config/phase5_motif_attribution.json \
  --attributions "$work_dir/dual_orientation_attributions.h5" --output-dir "$work_dir/output"
{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=$SLURM_JOB_ID"
  echo "node=$HOSTNAME"
  sha256sum containers/chrombpnet-upstream-20260914.sif "$source_h5" config/phase5_motif_attribution.json
} > "$work_dir/output/runtime.txt"
mv "$work_dir/output" "$output_dir"
rm "$work_dir/dual_orientation_attributions.h5"
rmdir "$work_dir"
(cd "$output_dir" && find . -maxdepth 1 -type f ! -name output_checksums.sha256 -print0 \
  | sort -z | xargs -0 sha256sum > output_checksums.sha256)
echo 'Motif/parent attribution overlap completed; scientific review remains pending.'
