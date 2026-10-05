# Parent-centered experimental planning draft

`parent_gain_control_groups.tsv.gz` has one row per ranked natural locus, with the original 500bp parent sequence and its proposed gain/loss-control sequences, coordinates, direct model effects, parent attribution, DA and donor summaries, repeat annotations and reference mappability. Sequence and coordinate hashes were checked against the integrated library; no experimental subset or oligo adapters were selected.

| Group composition | Loci | Unique inserts |
|---|---:|---:|
| Parent + gain + loss control | 938 | 2,814 |
| Parent + loss control | 8 | 16 |
| Parent only | 54 | 54 |
| Total | 1,000 | 2,884 |

All 938 gains have a matching natural parent and loss control. The 54 parent-only loci comprise 53 attribution-QC failures and one passing locus without a qualifying edit. Retaining those natural parents in the complete draft does not imply their sequences have approved optimization variants. The full detailed integrated library and all pairwise/leave-one-donor-out statistics remain in `../motif_attribution/ranked_library_with_parent_attribution.tsv.gz`.

The current 2,884 count describes unique inserts before adapters and barcodes. Reporter barcode multiplicity, biological replicates, additional assay controls and cloning design can change the number of synthesized constructs. The natural-locus specificity score is measured from ATAC; mutant specificity is not predicted. These are model-guided accessibility candidates, not experimentally confirmed enhancers.

Final synthesis preparation needs the assay/vector, construct budget, insert-length constraints, cloning/adaptor sequences and any prohibited sequence motifs, and barcode/replicate design. No new HPC training or motif job is scheduled while those inputs remain unspecified. A smaller experimental subset should retain complete parent/gain/control groups where possible; this file does not choose a subset automatically.

Reproduce with a new output directory:

```bash
.venv/bin/python scripts/prepare_phase5_experimental_groups.py \
  --output-dir reports/phase5/experimental_groups_new
```
