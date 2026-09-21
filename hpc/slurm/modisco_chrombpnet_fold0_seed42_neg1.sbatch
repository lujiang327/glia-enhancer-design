#!/usr/bin/env bash
#SBATCH --job-name=mg_cbmotif_f0_s42_n1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
cd "$project_root"
module load singularity

train_dir="models/chrombpnet/fold_0/seed_42_neg1"
attribution_dir="$train_dir/qc_attribution"
motif_dir="$attribution_dir/modisco"
profile_scores="$attribution_dir/scores/muller_chrombpnet_nobias.profile_scores.h5"
counts_scores="$attribution_dir/scores/muller_chrombpnet_nobias.counts_scores.h5"
motif_database="data/intermediate/vendor/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c/chrombpnet/data/motifs.meme.txt"

required=("$profile_scores" "$counts_scores" "$motif_database")
for path in "${required[@]}"; do
  [[ -s "$path" ]] || { echo "Missing full-model TF-MoDISco input: $path" >&2; exit 1; }
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
  scores="$attribution_dir/scores/muller_chrombpnet_nobias.${score_type}_scores.h5"
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
  "acceptance_rule": "Leading no-bias patterns should contain coherent biological TF-family motifs supported by Müller-glia biology; dominance by Tn5/enzyme-bias patterns, diffuse GC composition, or repeats requires investigation or rejection.",
  "scope": "Motif names alone do not establish exact TF identity. Review profile and count patterns, seqlet support, residual Tn5 signal, and paired-RNA support."
}
EOF

tar -czf "$motif_dir/chrombpnet_neg1_motif_reports.tar.gz" \
  -C "$motif_dir" profile_report counts_report
find "$motif_dir" -type f ! -name output_checksums.sha256 -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > "$motif_dir/logs/output_checksums.sha256"

echo "Full-model TF-MoDISco completed. Download reports for manual biological review."
