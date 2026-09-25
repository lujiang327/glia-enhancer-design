#!/usr/bin/env bash
#SBATCH --job-name=mg_phase3_consensus
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
output_dir="$scratch_root/consensus_peaks"
cd "$project_root"
mkdir -p logs/slurm
[[ ! -e "$output_dir" ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 1; }

python3 scripts/build_phase3_consensus_peaks.py \
  --peak-root "$scratch_root/peaks" \
  --cell-types config/phase3_cell_types.tsv \
  --chrom-sizes data/raw/reference/hg38.canonical.chrom.sizes \
  --blacklist data/raw/hg38-blacklist.v2.bed.gz \
  --output-dir "$output_dir" \
  --width 500 \
  --min-donors 2

find "$output_dir" -maxdepth 1 -type f ! -name output_checksums.sha256 -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > "$output_dir/output_checksums.sha256"
echo "Consensus peak universe built. Count-matrix QC remains required before differential testing."
