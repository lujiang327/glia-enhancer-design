# Candidate motif scan recovery and review

Job 63282882 scanned all 2,884 library sequences in 9m50s, but exited FAILED (1:0) during postprocessing. The parser discarded FIMO 4.11.2's commented legacy header (`#pattern name`, `sequence name`), then treated the first match as column names. The scan itself completed successfully before this exception. This was a format-handling bug, not a model failure.

The raw output was preserved in the cluster repository under ignored `data/intermediate/phase5_library_motifs_job_63282882/fimo_hits.tsv`. Its SHA256 is `834b2686c88caeff8335d17f7d5cef81899009ac52574138293349468bcb451a`. The cluster revision was `3d2ccb2f231bac7e0ec16342e283b059fc6359c5`; the pinned container SHA is `3eea58b1606fe1245ad70abe70ea6268cf10587ea88de52263c268d1952249d5`. Raw output was downloaded and verified, then postprocessed locally without rescanning. It and the full annotated site table remain in ignored local intermediate directories; compact review outputs are checked in.

The parser now accepts both legacy and modern headers and tests first-match preservation, empty output and invalid headers. Recovery requires an explicit raw-output hash. All 41 tests pass. The original failed Slurm accounting record remains failed; scientific output recovery is recorded separately, not rewritten as cluster-job success.

All 2,884 sequences and their original model, donor, DA, sequence and annotation columns were independently checked for preservation. All 1,154,982 reported matches obey the nominal p-value threshold, allowing printed rounding to 0.0001 (497 sites); none has a q-value. No duplicate site keys were found. Independent parent/variant site comparisons confirmed every gained/lost match overlaps its substituted base.

| Sequence role | Sequences | At least one gained site | At least one lost site | Any threshold crossing |
|---|---:|---:|---:|---:|
| Gain variant | 938 | 828 | 718 | 902 |
| Loss control | 946 | 634 | 921 | 937 |
| Natural parent | 1,000 | — | — | — |

The 36 gain variants and 9 loss controls without a threshold crossing still have direct model effects. A lack of threshold crossing does not establish a lack of motif-score change or motif involvement. Gains and losses can co-occur in the same variant. The median parent has 374 matches because this broad motif database contains redundant family models and enzyme motifs. Match counts are not independent regulatory sites or a new sequence ranking.

`library_motifs/ranked_library_with_motif_matches.tsv.gz` retains the ranked parent/variant library with sequence-match summaries. `edit_overlapping_motif_comparisons.tsv.gz` contains 59,828 edit-overlapping comparisons, including retained, gained and lost matches with strand, coordinates, parent/variant site scores and p-values where reported. A missing score means the site was below the reporting threshold, not that its score was zero. Direct model effects in both orientations accompany each comparison.

These are nominal FIMO sequence matches, not FDR-controlled binding evidence or causal attribution. RNA support and model-wide MoDISco motifs must not be assigned automatically to individual candidates. Per-site overlap with parent count/profile attribution is still pending; no variant attribution has been computed. Variant off-target specificity remains unavailable. The library remains a review draft, with adapter/barcode/cloning design and final experimental selection pending.

The next analysis should connect edit-overlapping matches to the existing count/profile attribution, examine representative gain/control pairs and ambiguous repeat-associated loci, and propose a manageable experimental subset. Retain natural parents and loss controls alongside proposed gains.

Reproduce local recovery (the output directory must be new):

```bash
.venv/bin/python scripts/scan_phase5_library_motifs.py \
  --config config/phase5_library_motifs.json \
  --output-dir data/intermediate/phase5_library_motifs_recovered \
  --reuse-fimo data/intermediate/phase5_library_motifs_recovery/fimo_hits.tsv \
  --reuse-fimo-sha256 834b2686c88caeff8335d17f7d5cef81899009ac52574138293349468bcb451a \
  --reuse-fimo-version 4.11.2
.venv/bin/python scripts/review_phase5_library_motifs.py \
  --results-dir data/intermediate/phase5_library_motifs_recovered \
  --output-dir reports/phase5/library_motifs
```
