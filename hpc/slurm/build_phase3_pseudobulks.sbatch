#!/usr/bin/env bash
#SBATCH --job-name=mg_phase3_pseudobulk
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --array=0-12%4
#SBATCH --output=logs/slurm/%x-%A_%a.out
#SBATCH --error=logs/slurm/%x-%A_%a.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
cd "$project_root"
mkdir -p logs/slurm "$scratch_root/pseudobulk"

row=$(awk -F'\t' -v line="$((SLURM_ARRAY_TASK_ID + 2))" 'NR == line {print; exit}' config/phase3_cell_types.tsv)
[[ -n "$row" ]] || { echo "Missing cell-type map row" >&2; exit 1; }
IFS=$'\t' read -r file_label published_label analysis_group priority <<< "$row"

fragments="$scratch_root/raw/GSM5866073_${file_label}_frags.tsv.gz"
assignments="$scratch_root/audit/$analysis_group/cell_assignments.tsv.gz"
output_dir="$scratch_root/pseudobulk/$analysis_group"
for path in "$fragments" "$assignments" \
  data/raw/reference/hg38.canonical.chrom.sizes \
  data/raw/hg38-blacklist.v2.bed.gz; do
  [[ -s "$path" ]] || { echo "Missing pseudobulk input: $path" >&2; exit 1; }
done
[[ ! -e "$output_dir" ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 1; }

python3 scripts/build_phase3_pseudobulks.py \
  --fragments "$fragments" \
  --assignments "$assignments" \
  --chrom-sizes data/raw/reference/hg38.canonical.chrom.sizes \
  --blacklist data/raw/hg38-blacklist.v2.bed.gz \
  --analysis-group "$analysis_group" \
  --published-label "$published_label" \
  --output-dir "$output_dir"

sha256sum "$output_dir/pseudobulk_summary.json" \
  > "$output_dir/summary.sha256"
echo "Completed Phase 3 pseudobulks for $analysis_group ($priority)."
