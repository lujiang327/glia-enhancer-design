#!/usr/bin/env bash
#SBATCH --job-name=retina_browser
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --array=0-12%3
#SBATCH --output=logs/slurm/%x-%A_%a.out
#SBATCH --error=logs/slurm/%x-%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Set scratch root}"
group=$(python3 -c 'import json,sys; print(json.load(open("config/browser_celltype_tracks.json"))["group_order"][int(sys.argv[1])])' "$SLURM_ARRAY_TASK_ID")
source_dir="$scratch_root/pseudobulk/$group"
output_dir="$scratch_root/browser_tracks_100bp/$group"
work_dir="data/intermediate/browser_celltype_job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
[[ -d "$source_dir" && ! -e "$output_dir" && ! -e "$work_dir" ]] || { echo 'Missing source or output already exists' >&2; exit 1; }
mkdir -p "$work_dir/input" "$scratch_root/browser_tracks_100bp"
cp "$source_dir/"*.fragments.bedpe.gz "$work_dir/input/"
module load singularity
extra_args=()
if [[ "$group" == MG ]]; then extra_args=(--mg-reference data/intermediate/training/muller.insertions.bw); fi
hpc/run_chrombpnet_container.sh python -u scripts/build_celltype_browser_tracks.py \
  --config config/browser_celltype_tracks.json --group "$group" --pseudobulk-dir "$work_dir/input" \
  --chrom-sizes data/raw/reference/hg38.canonical.chrom.sizes --output-dir "$work_dir/output" "${extra_args[@]}"
{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=$SLURM_JOB_ID"
  echo "node=$HOSTNAME"
  sha256sum config/browser_celltype_tracks.json containers/chrombpnet-upstream-20260914.sif
} > "$work_dir/output/runtime.txt"
mv "$work_dir/output" "$output_dir"
rm "$work_dir/input/"*.fragments.bedpe.gz
rmdir "$work_dir/input" "$work_dir"
(cd "$output_dir" && find . -maxdepth 1 -type f ! -name output_checksums.sha256 -print0 | sort -z | xargs -0 sha256sum > output_checksums.sha256)
echo "Completed observed browser tracks for $group; no off-target sequence model was trained."
