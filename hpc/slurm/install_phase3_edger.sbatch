#!/usr/bin/env bash
#SBATCH --job-name=mg_phase3_edger_install
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
runtime_root="$scratch_root/software/R-4.4.3-bioc-3.19-edger"
part_root="${runtime_root}.part.${SLURM_JOB_ID}"
cd "$project_root"
mkdir -p logs/slurm "$scratch_root/software"
[[ ! -e "$runtime_root" ]] || { echo "Refusing to overwrite: $runtime_root" >&2; exit 1; }
[[ ! -e "$part_root" ]] || { echo "Refusing to overwrite: $part_root" >&2; exit 1; }
mkdir -p "$part_root/library"

module load R/4.4.3
Rscript scripts/install_phase3_edger.R \
  --library "$part_root/library" \
  --manifest-dir "$part_root"

(cd "$part_root" && \
  find . -maxdepth 1 -type f ! -name runtime_checksums.sha256 -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > runtime_checksums.sha256)
printf 'runtime_ready\n' > "$part_root/.complete"
mv "$part_root" "$runtime_root"
echo "Phase 3 edgeR runtime installed at $runtime_root"
