#!/usr/bin/env bash
#SBATCH --job-name=mg_phase3_pooled_peaks
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH --array=0-12%3
#SBATCH --output=logs/slurm/%x-%A_%a.out
#SBATCH --error=logs/slurm/%x-%A_%a.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
cd "$project_root"
module load singularity
mkdir -p logs/slurm "$scratch_root/peaks"

row=$(awk -F'\t' -v line="$((SLURM_ARRAY_TASK_ID + 2))" 'NR == line {print; exit}' config/phase3_cell_types.tsv)
[[ -n "$row" ]] || { echo "Missing cell-type map row" >&2; exit 1; }
IFS=$'\t' read -r file_label published_label analysis_group priority <<< "$row"
output_dir="$scratch_root/peaks/$analysis_group/pooled"
[[ ! -e "$output_dir" ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 1; }

inputs=()
container_inputs=()
for donor in LGS1 LGS2 LGS3 LVG1; do
  input="$scratch_root/pseudobulk/$analysis_group/$donor.fragments.bedpe.gz"
  [[ -s "$input" ]] || { echo "Missing donor pseudobulk: $input" >&2; exit 1; }
  inputs+=("$input")
  container_inputs+=("/phase3/pseudobulk/$analysis_group/$donor.fragments.bedpe.gz")
done
mkdir -p "$output_dir"

PHASE3_SCRATCH_ROOT="$scratch_root" hpc/run_phase3_macs3.sh \
  macs3 callpeak \
  -t "${container_inputs[@]}" \
  -f BEDPE \
  -g hs \
  --keep-dup all \
  -q 0.01 \
  --call-summits \
  -n "${analysis_group}_pooled" \
  --outdir "/phase3/peaks/$analysis_group/pooled"

sha256sum "$output_dir/${analysis_group}_pooled_peaks.narrowPeak" \
  > "$output_dir/peak.sha256"
echo "Completed pooled peak call for $analysis_group ($priority)."
