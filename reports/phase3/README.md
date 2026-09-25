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

All 52 pseudobulks passed the depth and filtering review. See
`pseudobulk_qc_review.md`. The next permitted computation is reproducible peak
calling and common-universe count-matrix QC; differential testing is still
blocked.

## Reproducible peak-calling gate

The validated ChromBPNet image does not contain MACS3. Install MACS3 3.0.4 in
an isolated scratch environment, leaving the model runtime unchanged:

```bash
MACS_JOB=$(sbatch --parsable \
  --account=thahoang0 \
  --export=ALL,PHASE3_SCRATCH_ROOT="$PHASE3_SCRATCH_ROOT" \
  hpc/slurm/install_phase3_macs3.sbatch)
MACS_JOB=${MACS_JOB%%;*}
```

Submit donor and pooled calls after the installation succeeds, then build the
common universe only after both arrays finish:

```bash
DONOR_PEAK_JOB=$(sbatch --parsable \
  --account=thahoang0 \
  --dependency="afterok:${MACS_JOB}" \
  --export=ALL,PHASE3_SCRATCH_ROOT="$PHASE3_SCRATCH_ROOT" \
  hpc/slurm/call_phase3_donor_peaks.sbatch)
DONOR_PEAK_JOB=${DONOR_PEAK_JOB%%;*}

POOLED_PEAK_JOB=$(sbatch --parsable \
  --account=thahoang0 \
  --dependency="afterok:${MACS_JOB}" \
  --export=ALL,PHASE3_SCRATCH_ROOT="$PHASE3_SCRATCH_ROOT" \
  hpc/slurm/call_phase3_pooled_peaks.sbatch)
POOLED_PEAK_JOB=${POOLED_PEAK_JOB%%;*}

CONSENSUS_JOB=$(sbatch --parsable \
  --account=thahoang0 \
  --dependency="afterok:${DONOR_PEAK_JOB}:${POOLED_PEAK_JOB}" \
  --export=ALL,PHASE3_SCRATCH_ROOT="$PHASE3_SCRATCH_ROOT" \
  hpc/slurm/build_phase3_consensus_peaks.sbatch)
CONSENSUS_JOB=${CONSENSUS_JOB%%;*}

printf 'MACS_JOB=%s\nDONOR_PEAK_JOB=%s\nPOOLED_PEAK_JOB=%s\nCONSENSUS_JOB=%s\n' \
  "$MACS_JOB" "$DONOR_PEAK_JOB" "$POOLED_PEAK_JOB" "$CONSENSUS_JOB"
```

The universe retains a pooled cell-type peak only when at least two donors have
an overlapping peak. Reproducible summits are converted to nonoverlapping 500-bp
regions across all 13 cell types. The complete parameters are recorded in
`config/phase3_peak_calling.json`. Count-matrix and FRiP QC remain required before
differential accessibility.
