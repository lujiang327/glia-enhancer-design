#!/usr/bin/env bash
#SBATCH --job-name=mg_phase3_macs3
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
environment="$scratch_root/software/macs3-3.0.4"
cd "$project_root"
module load singularity
mkdir -p logs/slurm "$scratch_root/software"
[[ ! -e "$environment" ]] || { echo "Refusing to overwrite: $environment" >&2; exit 1; }

cleanup_failed_install() {
  if [[ -d "$environment" && ! -s "$environment/INSTALL_COMPLETE" ]]; then
    mv "$environment" "${environment}.failed-${SLURM_JOB_ID:-unknown}"
  fi
}
trap cleanup_failed_install ERR

singularity exec \
  --cleanenv \
  --no-mount hostfs \
  --bind "$project_root:/workspace" \
  --bind "$scratch_root:/phase3" \
  containers/chrombpnet-upstream-20260914.sif \
  bash -c '
    set -euo pipefail
    export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin
    python -m venv /phase3/software/macs3-3.0.4
    /phase3/software/macs3-3.0.4/bin/python -m pip install --no-cache-dir MACS3==3.0.4
    /phase3/software/macs3-3.0.4/bin/macs3 --version
    /phase3/software/macs3-3.0.4/bin/python -m pip freeze \
      > /phase3/software/macs3-3.0.4/requirements.freeze.txt
  '

completion_part="$environment/INSTALL_COMPLETE.part"
{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unset}"
  singularity --version
  sha256sum containers/chrombpnet-upstream-20260914.sif
  PHASE3_SCRATCH_ROOT="$scratch_root" hpc/run_phase3_macs3.sh macs3 --version
  sha256sum "$environment/requirements.freeze.txt"
} > "$completion_part"
mv "$completion_part" "$environment/INSTALL_COMPLETE"

echo "Installed isolated MACS3 3.0.4 environment."
