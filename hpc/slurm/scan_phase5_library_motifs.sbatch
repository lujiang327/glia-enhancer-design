#!/usr/bin/env bash
#SBATCH --job-name=mg_library_motifs
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from repository root}"
output_dir="${PHASE3_SCRATCH_ROOT:?Set scratch root}/phase5_library_motifs"
staging_dir="data/intermediate/phase5_library_motifs_job_${SLURM_JOB_ID}"
[[ ! -e "$output_dir" && ! -e "$staging_dir" ]] || { echo 'Output already exists' >&2; exit 1; }
module load singularity
hpc/run_chrombpnet_container.sh python -u scripts/scan_phase5_library_motifs.py \
  --config config/phase5_library_motifs.json --output-dir "$staging_dir"
{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=$SLURM_JOB_ID"
  echo "node=$HOSTNAME"
  sha256sum containers/chrombpnet-upstream-20260914.sif
} > "$staging_dir/runtime.txt"
mv "$staging_dir" "$output_dir"
(cd "$output_dir" && find . -maxdepth 1 -type f ! -name output_checksums.sha256 -print0 \
  | sort -z | xargs -0 sha256sum > output_checksums.sha256)
echo 'Library motif scan completed; scientific review pending.'
