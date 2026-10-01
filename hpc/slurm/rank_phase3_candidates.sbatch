#!/usr/bin/env bash
#SBATCH --job-name=mg_phase3_rank
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
results_dir="$scratch_root/differential_accessibility"
output_dir="$scratch_root/candidate_ranking"
matrix="$scratch_root/count_matrix/retina_consensus_peak_counts.tsv.gz"
normalization="reports/phase3/differential_accessibility_qc/normalization_factors.tsv"
refgene="data/raw/hg38.refGene.txt.gz"

cd "$project_root"
mkdir -p logs/slurm
for path in \
  "$results_dir/contrasts/MG_vs_rest.tsv.gz" \
  "$matrix" \
  "$normalization" \
  "$refgene" \
  config/phase3_candidate_ranking.json; do
  [[ -s "$path" ]] || { echo "Missing candidate-ranking input: $path" >&2; exit 1; }
done
[[ ! -e "$output_dir" ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 1; }

python3 scripts/rank_phase3_candidates.py \
  --results-dir "$results_dir" \
  --count-matrix "$matrix" \
  --normalization "$normalization" \
  --refgene "$refgene" \
  --config config/phase3_candidate_ranking.json \
  --output-dir "$output_dir"

(cd "$output_dir" && \
  find . -maxdepth 1 -type f ! -name output_checksums.sha256 -print0 \
    | sort -z | xargs -0 sha256sum > output_checksums.sha256)

echo "Candidate ranking completed. Review is required before ChromBPNet sequence scoring."
