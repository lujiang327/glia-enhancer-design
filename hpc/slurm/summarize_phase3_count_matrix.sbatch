#!/usr/bin/env bash
#SBATCH --job-name=mg_phase3_countqc
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
scratch_root="${PHASE3_SCRATCH_ROOT:?Submit with PHASE3_SCRATCH_ROOT set}"
cd "$project_root"
mkdir -p logs/slurm
module load singularity

image="containers/chrombpnet-upstream-20260914.sif"
[[ -s "$image" ]] || { echo "Missing container: $image" >&2; exit 1; }
[[ ! -e "$scratch_root/count_matrix" ]] || { echo "Refusing to overwrite count matrix" >&2; exit 1; }
[[ ! -e reports/phase3/count_matrix_qc ]] || { echo "Refusing to overwrite count-matrix QC report" >&2; exit 1; }

singularity exec \
  --cleanenv \
  --no-mount hostfs \
  --bind "$PWD:/workspace" \
  --bind "$scratch_root:/phase3" \
  "$image" \
  bash -c '
    export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin
    export MPLCONFIGDIR=/tmp/matplotlib-${SLURM_JOB_ID:-phase3}
    cd /workspace
    python scripts/summarize_phase3_count_matrix.py \
      --count-root /phase3/counts \
      --peaks /phase3/consensus_peaks/retina_consensus_peaks.bed \
      --cell-types config/phase3_cell_types.tsv \
      --pseudobulk-qc reports/phase3/pseudobulk_qc.tsv \
      --matrix-output-dir /phase3/count_matrix \
      --report-output-dir reports/phase3/count_matrix_qc
  '

echo "Count-matrix QC assembled. Differential accessibility remains blocked pending review."
