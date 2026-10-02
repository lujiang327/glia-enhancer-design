#!/usr/bin/env bash
#SBATCH --job-name=mg_locus_annot
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail
project_root="${SLURM_SUBMIT_DIR:?Submit from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
output_dir="$scratch_root/phase4_locus_annotation"
staging_dir="data/intermediate/phase4_locus_annotation_job_${SLURM_JOB_ID}"
cd "$project_root"
module load singularity
[[ ! -e "$output_dir" ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 1; }
[[ ! -e "$staging_dir" ]] || { echo "Refusing to overwrite: $staging_dir" >&2; exit 1; }
for path in reports/phase4/parent_scoring/parent_predictions.tsv.gz \
  reports/phase4/perturbation_full/best_passing_single_edits.tsv.gz \
  data/raw/reference/hg38.canonical.fa.fai config/phase4_locus_annotation.json; do
  [[ -s "$path" ]] || { echo "Missing annotation input: $path" >&2; exit 1; }
done
mkdir -p "$output_dir.runtime"
{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unset}"
  echo "node=$HOSTNAME"
  singularity --version
  sha256sum containers/chrombpnet-upstream-20260914.sif config/phase4_locus_annotation.json
} > "$output_dir.runtime/runtime.txt"

# No GPU is needed. The source snapshot is kept under ignored raw data.
hpc/run_chrombpnet_container.sh python -c 'import pandas, pyBigWig; print("Annotation dependencies available", flush=True)'
python3 -u scripts/download_phase4_annotation_assets.py --config config/phase4_locus_annotation.json
hpc/run_chrombpnet_container.sh \
  python -u scripts/annotate_phase4_loci.py \
  --config /workspace/config/phase4_locus_annotation.json \
  --output-dir "/workspace/$staging_dir"

mv "$staging_dir" "$output_dir"
mv "$output_dir.runtime/runtime.txt" "$output_dir/runtime.txt"
rmdir "$output_dir.runtime"
(cd "$output_dir" && \
  find . -maxdepth 1 -type f ! -name output_checksums.sha256 -print0 \
    | sort -z | xargs -0 sha256sum > output_checksums.sha256)
echo "Locus annotations completed. Review repeat and mappability context before library selection."
