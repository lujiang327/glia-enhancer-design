# Local training-input preparation

Run from the existing repository root. Complete the donor-resolved QC workflow
in the README first. These steps use the Mac CPU; they do not train a model or
perform differential accessibility. The private `config/greatlakes.json` stays
outside Git and is copied manually.

```bash
.venv/bin/python -m pip install -r config/preparation_requirements.lock.txt
.venv/bin/python scripts/download.py config/reference_assets.json
.venv/bin/python scripts/run_preparation.py reference > logs/prepare_reference.log 2>&1
.venv/bin/python scripts/run_preparation.py shift > logs/check_fragment_shift.log 2>&1
.venv/bin/python scripts/run_preparation.py tracks > logs/build_training_tracks.log 2>&1
.venv/bin/python scripts/run_preparation.py peaks > logs/training_peaks.log 2>&1
.venv/bin/python scripts/run_preparation.py regions > logs/training_regions.log 2>&1
.venv/bin/python scripts/run_preparation.py validate > logs/training_validation.log 2>&1
.venv/bin/python scripts/report_training_preparation.py > logs/preparation_report.log 2>&1
```

Execute in order and stop on any error. Reruns overwrite stage outputs; preserve
the run directory first when comparing versions. Inspect the reports under
`reports/preparation/` before using the inputs. Input validation is not a model
reliability assessment. `config/model_validation.json` must remain NOT_TRAINED
until actual model fitting and scientific evaluation are complete.

The saved `gc_review.json` records the current pilot acceptance judgment. Renew
that review after changing parameters or input data; report generation does not
automatically make a new scientific acceptance decision.

## Reference and signal conventions

The UCSC hg38 analysis-set FASTA is checked against the publisher MD5 and a saved
SHA256, then restricted to chr1–22 and X. Uppercasing preserves coordinate lengths
and hard-masked N bases. Mitochondrial, Y and alternate contigs are excluded
consistently from these model inputs and chromosome folds.

Each donor contributes one molecule per retained fragment record, regardless of
the duplicate-support column in the original source. Molecules in different cells
or donors with identical coordinates are retained. Donor fragments are never
deduplicated together. Concatenated eye samples can revisit a chromosome; sparse
counts are accumulated across those blocks before writing each donor bigWig.

Shift checks use the ATAC reference PWMs and pure shift calculations extracted
from ChromBPNet commit `eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c` (v1.0.1),
with the corresponding Modisco IC calculation. Sources, licenses and hashes are
saved under `scripts/vendor/`. Two independent 50,000-fragment samples per donor
are drawn deterministically. The detected native +4/−5 convention places both
one-base insertion records at the original start and end under the upstream
target +4/−4 convention (including its terminal BED minus-one conversion).
The earlier QC TSS proxy used a different endpoint convention and remains
explicitly labeled a proxy.

`muller.original_fragments.tsv.gz` retains native three-column fragment endpoints
for upstream fragment-input routines that detect shifts themselves. The derived
insertion BED and bigWigs are already corrected signals. Record the precise
chosen upstream training entry point before configuring the HPC job.

## Regions and background sampling

Genome-wide relaxed MACS3 peaks use insertion BED, `--keep-dup all --nomodel
--shift -75 --extsize 150 -p 0.01 --call-summits`. This uses MACS3 3.0.4 in place
of upstream MACS2; the earlier diagnostic QC peak set is not the training set.
The eventual approximately 1,000 candidates have no role in defining this set.

Inputs are 2,114 bp with 1,000-bp outputs and up to 500 bp jitter. Full jitter
contexts must fit the chromosome, contain ACGT only, and avoid the blacklist.
Background contexts additionally exclude all raw relaxed peak spans. Candidate
backgrounds are sampled on a 100-bp grid and matched within each chromosome
partition to GC rounded to 0.01: two negatives per training/validation peak and
one per test peak, without replacement. When exact buckets are exhausted, the
nearest available bucket is used. This bounded fallback differs from upstream
random-walk matching; every deviation is recorded and must be reviewed.

The initial 1,000-bp grid produced substantial high-GC mismatch and was rejected.
Its parameters, metrics and log are preserved under
`reports/preparation/initial_grid_1000/`; its large candidate and matching files
are archived under `data/intermediate/training/initial_grid_1000/`. A denser grid
increases candidate availability but allows neighboring backgrounds to overlap
within a split. It does not create additional independent genomic loci.

Five saved chromosome folds keep train, validation and test chromosomes disjoint.
Validation checks exact sequence and reverse-complement duplicates across splits;
it does not exhaustively rule out near-homology or shifted homologous windows.
The same biological donors contribute to chromosome splits, so donor-held-out
evaluation is still needed to assess donor generalization.

## Transfer after validation

Only after `reports/preparation/input_validation.json` reports PASS and GC
matching has been reviewed, transfer paths listed in
`config/training_transfer_files.txt`. Code, public configuration and small reports
travel through the existing Git repository. Large data and checkpoints stay out
of Git; the virtual environment is recreated separately for each platform.

Set your destination privately in your shell and run from the Mac repository:

```bash
rsync -avh --progress --files-from=config/training_transfer_files.txt ./ "$TRAINING_DESTINATION/"
```

`TRAINING_DESTINATION` must be your SSH host plus the remote project directory
(`user@host:/path/to/project`). Pull the same Git commit in that directory and
verify the transferred bytes there with:

```bash
python3 scripts/verify_training_transfer.py
```

Then configure a pinned Linux GPU environment and perform a short Great Lakes
bias/model smoke test before scheduling full fitting. Record the code commit,
container digest, CUDA/TensorFlow versions, parameters, seeds and Slurm logs.
An input-integrity PASS does not establish model reliability. Differential
accessibility and design remain blocked until scientific model validation passes.
