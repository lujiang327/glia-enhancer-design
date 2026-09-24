#!/usr/bin/env bash
#SBATCH --job-name=mg_phase3_audit
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --array=0-12%4
#SBATCH --output=logs/slurm/%x-%A_%a.out
#SBATCH --error=logs/slurm/%x-%A_%a.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
cd "$project_root"
mkdir -p logs/slurm "$scratch_root/audit"

row=$(awk -F'\t' -v line="$((SLURM_ARRAY_TASK_ID + 2))" 'NR == line {print; exit}' config/phase3_cell_types.tsv)
[[ -n "$row" ]] || { echo "Missing cell-type map row" >&2; exit 1; }
IFS=$'\t' read -r file_label published_label analysis_group priority <<< "$row"
fragments="$scratch_root/raw/GSM5866073_${file_label}_frags.tsv.gz"
output_dir="$scratch_root/audit/$analysis_group"
[[ -s "$fragments" ]] || { echo "Missing fragment file: $fragments" >&2; exit 1; }
[[ ! -e "$output_dir" ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 1; }

python3 scripts/audit_phase3_celltype_fragments.py \
  --fragments "$fragments" \
  --rna-barcodes data/raw/GSM5866081_barcodes.tsv.gz \
  --donor-metadata reports/qc/donor_metadata.tsv \
  --file-label "$file_label" \
  --published-label "$published_label" \
  --analysis-group "$analysis_group" \
  --output-dir "$output_dir"

sha256sum "$output_dir/audit_summary.json" "$output_dir/cell_assignments.tsv.gz" \
  > "$output_dir/output_checksums.sha256"
echo "Completed Phase 3 barcode audit for $analysis_group ($priority)."
