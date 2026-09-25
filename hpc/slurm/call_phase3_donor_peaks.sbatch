#!/usr/bin/env bash
#SBATCH --job-name=mg_phase3_donor_peaks
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-51%4
#SBATCH --output=logs/slurm/%x-%A_%a.out
#SBATCH --error=logs/slurm/%x-%A_%a.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
cd "$project_root"
module load singularity
mkdir -p logs/slurm "$scratch_root/peaks"

donors=(LGS1 LGS2 LGS3 LVG1)
group_index=$((SLURM_ARRAY_TASK_ID / 4))
donor_index=$((SLURM_ARRAY_TASK_ID % 4))
row=$(awk -F'\t' -v line="$((group_index + 2))" 'NR == line {print; exit}' config/phase3_cell_types.tsv)
[[ -n "$row" ]] || { echo "Missing cell-type map row" >&2; exit 1; }
IFS=$'\t' read -r file_label published_label analysis_group priority <<< "$row"
donor="${donors[$donor_index]}"
input="$scratch_root/pseudobulk/$analysis_group/$donor.fragments.bedpe.gz"
output_dir="$scratch_root/peaks/$analysis_group/donor_$donor"
[[ -s "$input" ]] || { echo "Missing donor pseudobulk: $input" >&2; exit 1; }
[[ ! -e "$output_dir" ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 1; }
mkdir -p "$output_dir"

PHASE3_SCRATCH_ROOT="$scratch_root" hpc/run_phase3_macs3.sh \
  macs3 callpeak \
  -t "/phase3/pseudobulk/$analysis_group/$donor.fragments.bedpe.gz" \
  -f BEDPE \
  -g hs \
  --keep-dup all \
  -q 0.01 \
  --call-summits \
  -n "${analysis_group}_${donor}" \
  --outdir "/phase3/peaks/$analysis_group/donor_$donor"

sha256sum "$output_dir/${analysis_group}_${donor}_peaks.narrowPeak" \
  > "$output_dir/peak.sha256"
echo "Completed donor peak call for $analysis_group $donor ($priority)."
