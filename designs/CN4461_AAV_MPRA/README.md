# CN4461 pilot cloning verification

Updated November 10, 2026. The exact published enhancer replacement is verified;
the user selected parents 140, 25 and 779 as the initial cloning pilot. All three
pass computational sequence checks. Some named feature boundaries remain
inferred. No orders have been placed.

Start with [the verification report](verification/CN4461_verification_report.md).

Current deliverables are separate:

- [Selected three-construct order package](pilot_order/README.md): three unchanged
  550-bp enhancer fragments; vendor workbook, CSV, FASTA and expected plasmids.
- [Pooled-barcode blueprint](pooled_barcode_design/POOLED_BARCODE_PLAN.md): proposed
  ten 16-bp tags per enhancer, screened codebook, cassette replacement strategy,
  engineering examples and required experimental QC. These are not initial-order
  items or a full-library synthesis release.

- `reference/`: unchanged official CN4461 and CN3044 files, published metadata,
  primary WPRE3 deposits and reference checksums.
- `verification/`: feature evidence, alignment, common arms, annotated reference,
  cassette maps, parameters and verification summary.
- `verification/pilot_review_only/`: three illustrative 500-bp natural parents,
  550-bp synthesis fragments, 4,738-bp simulated plasmids, review manifest and QC.
  The selection was subsequently accepted; these older files retain their
  historical review labels. Use `pilot_order/` for the current selected package.
- `preflight/`: earlier inventory of the larger parent/variant library;
  no full-library cloning design is released.

Reproduce from the repository root (Biopython 1.85 and Matplotlib):

```bash
MPLCONFIGDIR=/tmp/cn4461-matplotlib .venv/bin/python scripts/verify_cn4461_pilot.py
```

The script checks native sequence concordance, published enhancer identity,
actual enzyme cuts, directed overlap assembly, reporter translation, all native
feature sequences and every base outside the enhancer. It does not establish
physical assembly, ITR folding, packaging or experimental reporter activity.
Identical preserved BarcodeS cannot distinguish parents in a pooled barcode assay.
