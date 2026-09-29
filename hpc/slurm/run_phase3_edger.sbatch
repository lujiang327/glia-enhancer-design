#!/usr/bin/env bash
#SBATCH --job-name=mg_phase3_edger
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
runtime_root="$scratch_root/software/R-4.4.3-bioc-3.19-edger"
results_dir="$scratch_root/differential_accessibility"
report_dir="reports/phase3/differential_accessibility_qc"
matrix="$scratch_root/count_matrix/retina_consensus_peak_counts.tsv.gz"
cd "$project_root"
mkdir -p logs/slurm
for path in "$runtime_root/.complete" "$matrix" config/phase3_cell_types.tsv; do
  [[ -s "$path" ]] || { echo "Missing differential-analysis input: $path" >&2; exit 1; }
done
[[ ! -e "$results_dir" ]] || { echo "Refusing to overwrite: $results_dir" >&2; exit 1; }
[[ ! -e "$report_dir" ]] || { echo "Refusing to overwrite: $report_dir" >&2; exit 1; }

module load R/4.4.3
export R_LIBS_USER="$runtime_root/library"
Rscript scripts/run_phase3_edger.R \
  --count-matrix "$matrix" \
  --cell-types config/phase3_cell_types.tsv \
  --results-dir "$results_dir" \
  --report-dir "$report_dir"

(cd "$results_dir" && \
  find . -type f ! -name output_checksums.sha256 -print0 \
    | sort -z | xargs -0 sha256sum > output_checksums.sha256)
(cd "$report_dir" && \
  find . -maxdepth 1 -type f ! -name output_checksums.sha256 -print0 \
    | sort -z | xargs -0 sha256sum > output_checksums.sha256)
echo "Differential-model diagnostics completed. Candidate selection remains blocked pending review."
