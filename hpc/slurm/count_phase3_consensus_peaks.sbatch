#!/usr/bin/env bash
#SBATCH --job-name=mg_phase3_count
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --array=0-51%8
#SBATCH --output=logs/slurm/%x-%A_%a.out
#SBATCH --error=logs/slurm/%x-%A_%a.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
cd "$project_root"
mkdir -p logs/slurm "$scratch_root/counts"

group_index=$((SLURM_ARRAY_TASK_ID / 4 + 2))
donor_index=$((SLURM_ARRAY_TASK_ID % 4))
row=$(awk -F'\t' -v line="$group_index" 'NR == line {print; exit}' config/phase3_cell_types.tsv)
[[ -n "$row" ]] || { echo "Missing cell-type map row" >&2; exit 1; }
IFS=$'\t' read -r file_label published_label analysis_group priority <<< "$row"
donors=(LGS1 LGS2 LGS3 LVG1)
donor="${donors[$donor_index]}"
sample_id="${analysis_group}__${donor}"

fragments="$scratch_root/pseudobulk/$analysis_group/$donor.fragments.bedpe.gz"
peaks="$scratch_root/consensus_peaks/retina_consensus_peaks.bed"
output_dir="$scratch_root/counts/$sample_id"
for path in "$fragments" "$peaks"; do
  [[ -s "$path" ]] || { echo "Missing count input: $path" >&2; exit 1; }
done
[[ ! -e "$output_dir" ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 1; }

python3 scripts/count_phase3_consensus_peaks.py \
  --fragments "$fragments" \
  --peaks "$peaks" \
  --sample-id "$sample_id" \
  --analysis-group "$analysis_group" \
  --donor "$donor" \
  --output-dir "$output_dir"

echo "Completed consensus-peak counting for $sample_id ($published_label; $priority)."
