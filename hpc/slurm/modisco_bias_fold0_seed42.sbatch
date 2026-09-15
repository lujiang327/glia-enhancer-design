#!/usr/bin/env bash
#SBATCH --job-name=mg_biasmotif_f0_s42
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

# Supply the private account only at submission.
project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
cd "$project_root"
module load singularity

run_name="fold_0/threshold_0.5"
attribution_dir="models/bias/$run_name/seed_42/qc_attribution"
motif_dir="$attribution_dir/modisco"
profile_scores="$attribution_dir/scores/muller_bias.profile_scores.h5"
counts_scores="$attribution_dir/scores/muller_bias.counts_scores.h5"
motif_database="data/intermediate/vendor/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c/chrombpnet/data/motifs.meme.txt"

required=("$profile_scores" "$counts_scores" "$motif_database")
for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "Missing TF-MoDISco input: $path" >&2; exit 1; }
done
if [[ -e "$motif_dir" ]]; then
  echo "Refusing to overwrite existing output: $motif_dir" >&2
  exit 1
fi
mkdir -p "$motif_dir/logs" logs/slurm

{
  date --iso-8601=seconds
  echo "git_commit=$(git rev-parse HEAD)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unset}"
  echo "node=$HOSTNAME"
  singularity --version
  sha256sum containers/chrombpnet-upstream-20260914.sif
  sha256sum "${required[@]}"
} > "$motif_dir/logs/runtime.txt"

for score_type in profile counts; do
  scores="$attribution_dir/scores/muller_bias.${score_type}_scores.h5"
  hpc/run_chrombpnet_container.sh \
    env "OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}" \
    modisco motifs \
    -i "/workspace/$scores" \
    -n 50000 \
    -o "/workspace/$motif_dir/${score_type}_modisco.h5" \
    -w 500

  hpc/run_chrombpnet_container.sh \
    env "OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}" \
    modisco report \
    -i "/workspace/$motif_dir/${score_type}_modisco.h5" \
    -o "/workspace/$motif_dir/${score_type}_report" \
    -m "/workspace/$motif_database"
done

cat > "$motif_dir/logs/review_status.json" <<'EOF'
{
  "status": "MOTIF_REPORTS_READY_MANUAL_REVIEW_REQUIRED",
  "acceptance_rule": "Top profile motifs should resemble Tn5/enzyme bias or unstructured repeats; obvious transcription-factor motifs cause rejection or further investigation.",
  "scope": "No automated biological acceptance decision is made from motif names alone."
}
EOF

tar -czf "$motif_dir/bias_motif_reports.tar.gz" \
  -C "$motif_dir" profile_report counts_report
find "$motif_dir" -type f ! -name output_checksums.sha256 -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > "$motif_dir/logs/output_checksums.sha256"

echo "TF-MoDISco completed. Download the report archive for manual motif-leakage review."
