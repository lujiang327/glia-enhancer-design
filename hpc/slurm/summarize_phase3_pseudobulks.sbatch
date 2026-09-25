#!/usr/bin/env bash
#SBATCH --job-name=mg_phase3_pb_summary
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
cd "$project_root"
mkdir -p logs/slurm reports/phase3
output_json="reports/phase3/pseudobulk_qc.json"
output_tsv="reports/phase3/pseudobulk_qc.tsv"
[[ ! -e "$output_json" && ! -e "$output_tsv" ]] || {
  echo "Refusing to overwrite Phase 3 pseudobulk QC" >&2
  exit 1
}

python3 scripts/summarize_phase3_pseudobulk_qc.py \
  --pseudobulk-root "$scratch_root/pseudobulk" \
  --output-json "$output_json" \
  --output-tsv "$output_tsv"

sha256sum "$output_json" "$output_tsv" \
  > reports/phase3/pseudobulk_qc.sha256
echo "Phase 3 pseudobulk QC completed. Review required before peak calling."
