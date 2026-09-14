# Glia enhancer design — model feasibility first

Start with [the measured Müller glia QC report](reports/phase1_qc_report.html)
or its [Markdown version](reports/phase1_qc_report.md).
The pipeline does **not** run differential accessibility or sequence design before
ChromBPNet validation. The [workflow specification](docs_workflow.md) documents the
full project and its scientific validation gate.

## Current implementation

- Public data inventory and checksum-verified selective download.
- Full Müller glia fragment audit, barcode-collision detection and explicit exclusions.
- Published-label RNA marker check without manual relabeling.
- Chromosome/blacklist filtering, pooled and donor TSS proxies, fragment lengths.
- Four identifiable donor pseudobulks, broad-bin donor consistency, diagnostic MACS3
  accessible peaks and FRiP. These are QC peak calls, not differential tests.
- Figures, measured tables, logs, pinned QC dependencies and five chromosome folds.
- Verified reference and donor/pooled insertion tracks, 256,248 genome-wide
  training peaks and GC-matched backgrounds; all five input-integrity checks pass.
- A fail-closed model-validation gate; later analyses are not implemented or run.

No trained ChromBPNet checkpoint exists yet. Local training-input preparation is
documented in [the preparation guide](docs/training_preparation.md). The GPU
runtime, bias model and ChromBPNet fitting/validation remain separate next steps.
See [the input-preparation report and GC plots](reports/preparation/README.md).

## Mac to Great Lakes

Use this existing repository, `https://github.com/lujiang327/glia-enhancer-design.git`,
for code, configuration, documentation and reviewable QC results. Prepare changes
on the Mac, push them to this repository, then pull the same commit on Great Lakes.
Record the commit used for each training run.

Private connection, account and directory settings are saved locally in
`config/greatlakes.json`. This file is excluded from Git and must be synchronized
manually between machines. Copy [the template](config/greatlakes.example.json)
to that path when configuring a new checkout. Singularity 4.4.1 was
confirmed by the user after `module load singularity`; the GPU training environment
has not yet been validated. These settings do not submit jobs or provision resources.

Git does not carry the ignored raw data, donor pseudobulks, local virtual environment,
container images or checkpoints. Recreate/download large inputs from the manifests
or transfer them separately, preserving their paths and verifying checksums. A
clone alone is therefore not sufficient to start training. After input validation,
the preparation guide explains how to transfer and verify the large inputs.

## Reproduce QC

Run from the repository root. Python 3.9.6 was used locally. The selective source
manifest contains sizes and SHA256 values; the four primary data assets also
match independently published HCA SHA256 values. Public endpoints may require
network-enabled execution. Upstream author code is retained for provenance, not
executed. Approximately 6 GB of local storage is used by this QC run.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r config/preparation_requirements.lock.txt
python3 scripts/download.py config/assets.json
mkdir -p reports/qc/peaks logs
python3 scripts/audit_fragments.py --fragments data/raw/GSM5866073_Mullerglia_frags.tsv.gz --rna-barcodes data/raw/GSM5866081_barcodes.tsv.gz --out reports/qc > logs/fragment_audit.log 2>&1
python3 scripts/marker_qc.py > logs/marker_qc.log 2>&1
python3 scripts/pool_qc.py --fragments data/raw/GSM5866073_Mullerglia_frags.tsv.gz --blacklist data/raw/hg38-blacklist.v2.bed.gz --refgene data/raw/hg38.refGene.txt.gz --sizes data/raw/hg38.chrom.sizes --out reports/qc > logs/pool_qc.log 2>&1
python3 scripts/donor_qc.py > logs/donor_qc.log 2>&1
.venv/bin/macs3 callpeak -t reports/qc/muller.filtered.bedpe -f BEDPE -g hs --keep-dup all -q 0.01 -n muller_diagnostic --outdir reports/qc/peaks > logs/macs3_qc.log 2>&1
python3 scripts/peak_qc.py > logs/peak_qc.log 2>&1
.venv/bin/python scripts/plot_qc.py > logs/plot_qc.log 2>&1
python3 scripts/report_qc.py > logs/report_qc.log 2>&1
.venv/bin/python -m unittest discover -s tests
python3 scripts/verify_outputs.py
```

Commands are sequential because each completed output is a prerequisite for the
next step. Large outputs are overwritten on a rerun: copy or version the run
directory first if preserving multiple runs. All steps fail on malformed inputs
or missing prerequisites; a partial download is never accepted as complete.

```bash
python3 scripts/model_gate.py
```

This currently exits with status 2, intentionally blocking phases 3–5. Passing
requires a validated checkpoint, matching checksum, all scientific checks marked
pass, and existing validation evidence. The record does not itself perform model
evaluation; it captures the outcome of the specified validation workflow.

## Saved artifacts

| Location | Contents |
|---|---|
| `data/raw/` | Unmodified fragment/RNA downloads and hg38 QC references |
| `data/metadata/` | GEO, HCA, CELLxGENE snapshots, complete file inventories, archived author code and model documentation |
| `data/intermediate/donor_pseudobulks/` | One filtered three-column MACS BEDPE per donor; no invented cell IDs |
| `reports/qc/` | Barcode audit, unambiguous cell mapping, marker summaries, depth/TSS/FRiP tables, QC bin counts and peaks |
| `reports/figures/` | PNG and SVG figures from measured data |
| `config/` | Verified source manifest, donor mapping, package lock, training plan, folds, QC recommendation and validation gate |
| `logs/` | Dependency install, data verification, QC, peak calling, figures and tests |

Full donor counts remain unknown for ambiguous barcodes. The report always
distinguishes the complete published pool from the usable donor-resolved subset.
Raw/large intermediate files are deliberately excluded from Git but retained
locally. Do not confuse future genomic regulatory classes with proven enhancers.
