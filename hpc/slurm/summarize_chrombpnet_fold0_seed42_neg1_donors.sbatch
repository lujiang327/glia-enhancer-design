#!/usr/bin/env bash
#SBATCH --job-name=mg_cbdonor_summary_f0_s42_n1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/slurm/%x-%j.out
#SBATCH --error=logs/slurm/%x-%j.err

set -euo pipefail

project_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
cd "$project_root"
module load singularity

qc_root="models/chrombpnet/fold_0/seed_42_neg1/qc_donor_sensitivity"
summary_json="$qc_root/donor_qc_summary.json"
summary_tsv="$qc_root/donor_qc_metrics.tsv"
if [[ -e "$summary_json" || -e "$summary_tsv" ]]; then
  echo "Refusing to overwrite an existing donor summary" >&2
  exit 1
fi
for donor in LGS1 LGS2 LGS3 LVG1; do
  for label in chrombpnet scaled_bias; do
    path="$qc_root/$donor/evaluation/muller_${label}_metrics.json"
    [[ -s "$path" ]] || { echo "Missing donor metric input: $path" >&2; exit 1; }
  done
done

hpc/run_chrombpnet_container.sh \
  python scripts/summarize_chrombpnet_donor_qc.py \
  --qc-root "/workspace/$qc_root" \
  --donor-summary /workspace/reports/qc/donor_summary.tsv \
  --output-json "/workspace/$summary_json" \
  --output-tsv "/workspace/$summary_tsv"

sha256sum "$summary_json" "$summary_tsv" > "$qc_root/donor_summary_checksums.sha256"
echo "Donor-sensitivity summary completed. Scientific review remains required."
