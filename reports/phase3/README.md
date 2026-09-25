# Phase 3 donor-by-cell-type audit

This directory contains the compact, versioned QC results used to decide whether
the published GSE196235 cell-type fragment pools can support donor-blocked
differential accessibility. Differential testing must not start until this audit
has been reviewed.

The 13 fragment pools are staged under `PHASE3_SCRATCH_ROOT/raw` on Great Lakes.
`audit_phase3_celltypes.sbatch` retains the published cell-type labels and assigns
a bare ATAC barcode to a donor only if that barcode matches exactly one prefixed
RNA barcode. Barcodes matching zero or multiple RNA samples are excluded from
donor-level pseudobulks. The audit reports unique fragment rows and the fifth-column
read-support sum separately, including canonical chr1-22/X totals.

Submit from the repository root:

```bash
export PHASE3_SCRATCH_ROOT=/scratch/thahoang_root/thahoang0/lujiang/chrombpnet_mg_phase3

AUDIT_JOB=$(sbatch --parsable \
  --account=thahoang0 \
  --export=ALL,PHASE3_SCRATCH_ROOT="$PHASE3_SCRATCH_ROOT" \
  hpc/slurm/audit_phase3_celltypes.sbatch)
AUDIT_JOB=${AUDIT_JOB%%;*}

SUMMARY_JOB=$(sbatch --parsable \
  --account=thahoang0 \
  --dependency="afterok:${AUDIT_JOB}" \
  --export=ALL,PHASE3_SCRATCH_ROOT="$PHASE3_SCRATCH_ROOT" \
  hpc/slurm/summarize_phase3_celltypes.sbatch)
SUMMARY_JOB=${SUMMARY_JOB%%;*}

printf 'AUDIT_JOB=%s\nSUMMARY_JOB=%s\n' "$AUDIT_JOB" "$SUMMARY_JOB"
```

Per-cell-type assignments and summaries remain in scratch under `audit/`. The
dependent summary job writes `donor_celltype_audit.json`,
`donor_celltype_audit.tsv`, and their checksum file here. Review zero-cell and
low-cell donor strata, donor balance, retained canonical read support, and mapping
loss before constructing donor-by-cell-type pseudobulks or calling consensus peaks.

The completed audit passed this gate with caveats. See
`donor_celltype_audit_review.md` for the decision and the required treatment of
low-cell strata and rod abundance during downstream analysis.

## Donor pseudobulk gate

After accepting the donor-by-cell-type audit, build the 52 filtered donor
pseudobulks and submit the dependent QC summary:

```bash
export PHASE3_SCRATCH_ROOT=/scratch/thahoang_root/thahoang0/lujiang/chrombpnet_mg_phase3

PSEUDOBULK_JOB=$(sbatch --parsable \
  --account=thahoang0 \
  --export=ALL,PHASE3_SCRATCH_ROOT="$PHASE3_SCRATCH_ROOT" \
  hpc/slurm/build_phase3_pseudobulks.sbatch)
PSEUDOBULK_JOB=${PSEUDOBULK_JOB%%;*}

PSEUDOBULK_SUMMARY_JOB=$(sbatch --parsable \
  --account=thahoang0 \
  --dependency="afterok:${PSEUDOBULK_JOB}" \
  --export=ALL,PHASE3_SCRATCH_ROOT="$PHASE3_SCRATCH_ROOT" \
  hpc/slurm/summarize_phase3_pseudobulks.sbatch)
PSEUDOBULK_SUMMARY_JOB=${PSEUDOBULK_SUMMARY_JOB%%;*}

printf 'PSEUDOBULK_JOB=%s\nPSEUDOBULK_SUMMARY_JOB=%s\n' \
  "$PSEUDOBULK_JOB" "$PSEUDOBULK_SUMMARY_JOB"
```

Each pseudobulk retains one count per fragment-file row after canonical-chromosome,
coordinate, and blacklist filtering. Fifth-column read support is recorded for
QC and is not used to duplicate fragments. Review `pseudobulk_qc.json` and
`pseudobulk_qc.tsv` before peak calling. Differential accessibility remains
blocked at this gate.
