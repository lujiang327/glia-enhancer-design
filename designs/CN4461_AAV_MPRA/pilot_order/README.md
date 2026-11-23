# Initial cloning pilot: three selected natural parents

November 10, 2026. The user selected the previously simulated parents 140,
25 and 779 for the initial cloning pilot. Their nucleotide sequences and
cloning arms are unchanged from the sequence-verified simulations.

For the synthesis vendor:

- **CN4461_three_construct_pilot_order.xlsx**: concise three-row order sheet
  plus a separate assembly/evidence review tab.
- **pilot_order.csv** and **pilot_order.fasta**: complete exact 550-bp fragments.
  These are alternative exports of the same three sequences, not extra items.
- **MG_P*_expected_plasmid.gbk**: expected final plasmids if the vendor performs
  cloning. Each is 4,738 bp and retains the original reference BarcodeS.

Request linear double-stranded DNA with no additional adapters or substitutions.
Ask the vendor to quote its standard minimum synthesis scale; yield, delivery
format, price and lead time have not been agreed. No purchase has been made.
Synthesis feasibility/sequence complexity and physical cloning success remain
vendor/laboratory checks.

If cloning is outsourced, agree how the vendor will receive a sequence-verified
physical CN4461 backbone; the digital reference file is not physical material.

## Exact arms, orientation and junctions

Both flanks below are written 5′→3′ in the supplied enhancer/SYFP2 direction:

```text
Left:  GTTCCTGCGGCCGCACGCGT            20 bp; CN4461 125–144
Right: GAGCTCGGGCTGGGCATAAAAGTCAGGGCA  30 bp; CN4461 727–756

Fragment = left + exact 500-bp parent + right = 550 bp
```

The right synthesis flank is appended as shown, not reverse-complemented.
The supplied parent is also used as shown, without reverse-complementing it.
Fragment bases 1–20 are the left flank, 21–520 the exact parent, 521–550 the
right flank. `pilot_sequence_audit.json` saves each junction and sequence hash.

MluI cuts the reference plus strand after 139 and the opposite strand after 143;
SacI cuts after 731 and727, respectively. The retained plus-strand overlap is
15 bp on the left (125–139) and 25 bp on the right (732–756). The remaining five
bases of each synthesis flank restore the complete restriction sites during
homology assembly. Estimated overlap Tm values are 54.5°C and 63.2°C under the
recorded assumptions. Complete arms are unique; the left arm shares a 14-bp
partial tract elsewhere and includes six terminal ITR bases. This is a physical
assembly consideration, not an observed cloning result.

Prepare the backbone separately with MluI/SacI and remove the original enhancer
fragment. Use two-fragment homology assembly; **do not digest the synthetic
inserts**. Parent 140 contains an internal SacI at enhancer 392–397, while parents 25
and 779 do not. This does not prevent the proposed separate-backbone digestion,
but changes any subsequent SacI diagnostic pattern.

## Verified details versus review items

Verified digitally: exact parent sequences; arm identity, lengths and direction;
restored junctions; intended 500-bp replacement; unchanged non-enhancer sequence;
native SYFP2 translation; 4,738-bp plasmid lengths; 2,141-bp long-annotation ITR
spans (2,130 bp with the alternative endpoint). Each exported sequence is checked
against the previous simulation and the authoritative reference.

Annotation review: exact minBG endpoints, WPRE3 terminal convention, original
BarcodeS name mapping, and overlapping right-ITR endpoints. The complete
surrounding sequences are retained despite those annotation limits.

Required physical QC: full insert and both junctions, entire plasmid identity,
reporter/barcode/polyA sequence, ITR integrity and absence of parental enhancer.
These are not yet measured. Vendor synthesis/cloning certificates do not establish
MG specificity, AAV genome integrity or reporter performance.

**All three retain AGCAATTA as the same existing barcode. Keep initial cloning
pilot samples separate; these tags cannot identify parents in a pooled readout.**
The separate `../pooled_barcode_design/` engineering files are not initial-order
items and must not be substituted into this pilot order.

Primary reference: https://www.addgene.org/236600/sequences/ (sequence472831).
Boundary evidence and confidence annotations: `../verification/CN4461_verification_report.md`.
Reproduce with `scripts/prepare_cn4461_pilot_and_barcodes.py`; the spreadsheet
builder is `scripts/create_cn4461_vendor_sheet.mjs` using Artifact Tool.
