# Pooled-screen barcode blueprint, separate from the cloning pilot

November 10, 2026. **Engineering design only.** Initial order items remain the
three unmodified-barcode 550-bp enhancer fragments in `../pilot_order/`.
This folder contains a candidate barcode codebook, replacement design, 30 digital
examples and QC plan. No eventual enhancer panel has been finalized or ordered.

## Barcode architecture and replication

Use **ten designed, distinct 16-bp transcribed barcodes per enhancer sequence**.
Allocate distinct barcodes to controls and every gain/loss variant as separate
sequence entities. Never reuse a tag across enhancer identities in the same
experiment. Ten is a proposed starting design, not a proven optimum: it allows
assessment of tag effects and some dropout. Aim to retain at least five usable
tags per enhancer in each biological sample; failure triggers review of power
and coverage, not silent acceptance of a single favorable tag.

`barcode_codebook_4000.csv/.fasta` supplies 4,000 reproducibly generated tags:
capacity for 400 enhancer/control sequences at ten tags each. For example,
300 enhancers plus 20 control sequences require 3,200 tags, leaving 800 reserve.
These are planning numbers, not a finalized selection or a synthesis release.
No full-library assignments are made. The first 30 tags are used only for
engineering examples on the three selected parents.

Design filters and actual exhaustive checks are in `codebook_QC.json`:

- Exactly 16 A/C/G/T bases; GC 7–9 bases (43.75%–56.25%).
- At most two identical bases in a row; no perfect six-base dinucleotide repeat.
- **Minimum pairwise Hamming distance 5**, including comparisons to reverse
  complements and each tag's own reverse complement.
- Avoid the defined cloning-site motifs and two common polyA hexamers within
  tags; check the actual upstream/downstream junctions for new sites as well.
- Record cycle-specific A/C/G/T frequencies; sequence with a validated read
  structure and appropriate library diversity rather than assuming the designed
  DNA alone guarantees sequencing performance.

These constraints do not prove regulatory neutrality, remove every possible
RNA-processing motif, or establish an indel-correcting code. The16-bp length is
a deliberate change from the published 8-bp tag and requires readout validation.
Preserving the same original tag in every construct is unsuitable for pooling.

**Decoding:** locate the exact known flanks, extract exactly 16 high-quality bases
in a predefined orientation, then use the expected sample whitelist. Accept an
exact tag or a unique match with at most one substitution; reject ambiguous,
wrong-length/indel, unexpected-flank or low-quality calls. Distance5 permits
stronger theoretical substitution correction, but the conservative decoder
does not use that as a reason to accept poor reads. A fixed tag is not a UMI:
use the assay's separate UMI for molecular counting where available.

Multiple tags are **technical observations**, not ten biological replicates.
Use biological replicates and donor/animal information for inference. Compare
RNA to matched vector-DNA abundance using library-size normalization and a
count model with DNA exposure/offset; assess tag effects and enhancer-level
agreement. Do not rank solely on raw RNA, one maximum tag or ratios from absent
DNA denominators. Report RNA and DNA counts, usable tags and dispersion.

## Replacing the existing BarcodeS sequence

The candidate existing tag is reference **1832–1839, AGCAATTA**, strongly supported
by CN3044 alignment and the published 8-bp barcode readout pattern. Its original
name-to-sequence mapping remains a depositor annotation question. The physical
sequence is known and can be targeted without inventing any bases.

```text
CN4461 original:  CTCGAG | AGCAATTA         | CGAGAGATCTTCGAC
Pooled design:    CTCGAG | NEW_16BP_BARCODE | CGAGAGATCTTCGAC
```

Retain every other nucleotide. Do not add a second barcode beside the old one,
delete the constant readout sequence or modify BGHpA to create room.

**Do not use XhoI as a unique linearization site:** reference XhoI starts are
205,399,593 and1826. Three sites are inside the original repeated enhancer.
Its use in the historical barcode method does not establish uniqueness for this
construct or our eventual candidate set.

The tested alternative is **EcoRV/PmlI replacement of a short cassette**:

| Detail, all reference coordinates 1-based inclusive | Verified design |
|---|---|
| EcoRV recognition |1574–1579, GAT/ATC; blunt cut after 1576 on both strands |
| PmlI recognition |2061–2066, CAC/GTG; blunt cut after 2063 on both strands |
| Unique in original and three selected parents |Yes, one site each; must rescreen every eventual sequence |
| Plus-strand sequence removed by cutting |1577–2063,487 nt |
| Cassette between complete restored sites |1580–2060,481 bp original; 489 bp after8→16 tag replacement |
| Left synthesis flank,5′→3′ |ACGGCGCGCCGCGGCCGCGAATTCGATATC,30 bp, 1550–1579 |
| Right synthesis flank,5′→3′ |CACGTGCGGACCGAGCGGCC,20 bp, 2061–2080 |
| Retained plus-strand overlaps |27 bp left;17 bp right |
| Synthetic cassette fragment | 30 + 489 + 20 = **539 bp** |

The literal cassette is `reference[1580–1831] + new 16 bp + reference[1840–2060]`.
Thus WPRE/linker, the entire known polyA sequence and its context are copied
unchanged from CN4461, not reconstructed from guessed feature definitions.
Both ITRs remain in the retained backbone; no whole-backbone PCR is required.
Digest only the backbone, then use homology assembly with the undigested
synthetic cassette and restore both complete sites. Screen cassette synthesis
complexity and actual assembly with the vendor/lab before scale-up.

The **30 examples** in `engineering_examples/` are three parents × ten tags.
All pass independent directed reconstruction from the actual cut backbone.
Every nonbarcode pilot base, original 500-bp parent and protected native feature
is unchanged. Expected plasmids are **4,746 bp**; long-annotation ITR span
**2,149 bp**, alternative short span2,138 bp. New tag coordinates in these
MG-parent plasmids are1750–1765; BGHpA1780–1891 and right ITR2009–2149.

`engineering_barcode_cassettes_NOT_ORDER.csv` is a design example, not a second
order sheet. `engineering_examples_expected_associations.tsv` is an **expected**
dictionary, not measured evidence of association. Raw low-error long reads or
validated consensus reads must confirm it experimentally.

For a few hundred enhancers, every new sequence must be screened for internal
EcoRV/PmlI sites, including parent/variant changes. A failed site screen requires
a separately verified local strategy or vendor targeted assembly; do not silently
digest through an enhancer or amplify the full ITR-bearing plasmid by default.
The30-example result does not establish compatibility for unselected sequences.

## Building and validating the enhancer–barcode dictionary

Maintain an immutable design table containing construct ID, full enhancer sequence
and hash, barcode ID/sequence, parent/variant/control relationship, expected
plasmid hash, clone/well, packaging unit and batch. Every later change creates a
new version; assignments are not recovered from abundance correlations.

1. **At plasmid construction:** for the first engineering panel, verify clones
   individually. For any larger pool, use reads spanning the **entire enhancer
   and barcode on the same molecule**, plus intervening cassette/junctions.
   Barcode-only short reads demonstrate abundance, not enhancer linkage.
2. **Build the measured dictionary:** align complete 500-bp candidate and tag
   flanks to the intended catalog; call orientation and structure, require one
   unambiguous enhancer/tag assignment and compare with the design dictionary.
   Record intended, wrong, unknown, mixed and incomplete molecules separately.
   Do not discard unexpected combinations before calculating discordance.
3. **Control preparation artifacts:** use PCR-free linkage sequencing where
   practical, matched plasmid controls and independently prepared controls.
   If amplification is required, quantify its contribution to chimera formation;
   do not equate PCR read count with independent original molecules. Near-identical
   parent/variant sequences need sufficient accuracy or consensus to distinguish
   the actual edited bases; ambiguous alignments remain unresolved.
4. **After AAV packaging:** repeat same-molecule enhancer→barcode linkage and
   genome-structure QC on protected packaged DNA, not just residual plasmid.
   Include a residual-plasmid contamination control and report complete versus
   truncated/structurally altered molecules. Long-read physical coverage,
   orientation and full-span callability must be reported before purity claims.
5. **In the assay:** barcode RNA reads use the measured dictionary. Sequence
   matched input and tissue/cell-type vector DNA, not only the plasmid pool;
   variable transduction and packaging are distinct from transcription.

Source evidence: pooled AAV packaging can break enhancer/tag linkage despite a
correct input plasmid dictionary. Additional barcodes do not repair that failure.
Published work also shows that individual packaging does not eliminate all
in-vivo multiplex noise, including possible transcriptional crosstalk.
See [Lalanne et al., Nature Biotechnology 2026](https://www.nature.com/articles/s41587-026-03097-1)
and [Hunker et al., Nature Communications 2026](https://www.nature.com/articles/s41467-026-72147-8).
Primary source XML is saved under `sources/` (the Lalanne XML is the earlier
preprint, explicitly labeled; the link points to its peer-reviewed publication).

## Packaging recommendation and staged QC

Start with **separate enhancer identities during packaging** and pool only
sequence/QC-qualified viral preparations. If the ten barcodes of one enhancer
are packaged together, swaps between those tags do not change enhancer identity,
but individual tag abundance/effects and genome structure still require QC.
An initial tag-effect comparison should use a small individually traceable subset.
Separate identities can require hundreds of production units: cost and production
capacity are an explicit scale-up constraint, not solved by this codebook.

If pooled packaging is needed for feasibility, benchmark it on a small barcoded
panel against separately packaged/mixed controls before committing the large
library. Include production replicates. Do not assume a dilution or capsid choice
eliminates swaps. No packaging choice here is a substitute for measuring linkage.

Before ordering the large pool, require these stages:

| Stage | Required evidence before advancing |
|---|---|
| Initial three-clone pilot |Complete insert/junctions; whole-plasmid identity; existing tag and both ITRs; no original enhancer contamination. |
| Small barcode engineering test |16-bp cassette replacement works; full sequence verified; tag readout and extraction work; tags do not materially alter expression, RNA processing or cell specificity. |
| Dictionary and representation |Every retained enhancer has multiple verified tags; expected→measured associations, counts, dropout and skew reported; no ambiguous reused tags. |
| Packaging benchmark |Protected-genome linkage/structure; residual-plasmid control; separately packaged/mixed comparator; packaging batch reproducibility. |
| Retinal functional test |MG and major retinal off-target readouts, vector-DNA normalization, biological replicates; compare individual and pooled behavior and delivery-dose effects. |
| Production release |Frozen catalog, assignments and reference files; explicit batch acceptance criteria; coverage/power and production budget; only qualifying identities carried forward. |

Controls to allocate fresh tags to: minimal-promoter-only, GC/length-matched
negative/shuffled sequences, appropriate motif-loss controls, and a supplied
experimentally validated MG-positive reporter. Do not label the current
striatal AiE0779m enhancer an established MG-positive control. Include no-template
and no-RT assay controls, cross-sample contamination controls, and single-vector
versus mixed-vector comparisons to reveal pooling/crosstalk effects. These
controls are proposed, not added to the initial three-fragment order.

**Proposed project acceptance targets, not established universal standards:**
all three physical pilot sequences must match their references; no verified
ITR/cassette alteration is accepted. For linkage, target at least 100 callable
full-span molecules per retained barcode and ≥1,000 per enhancer (400 identities
would require at least 400,000 callable molecules before overhead/dropout).
Require ≥99% observed intended linkage per enhancer and a 95% Wilson lower bound
≥98% in plasmid and packaged libraries; report both per-enhancer and per-tag
results and unknown/un-callable fractions. Insufficient coverage is inconclusive.
These stringent targets should be pre-agreed with the sequencing/packaging team;
failure means redesign/repackage/reduce scope, not accepting high-swapping batches.
Do not claim a passing global average rescues a failing enhancer.

Minimum five usable barcodes per enhancer per biological sample is a proposed
retention target, not proof of adequate power. Require a preplanned biological
replicate and coverage analysis, concordance across tags and batches, and direct
specificity evidence. Multiple tags, accurate linkage and successful cloning do
not themselves demonstrate an MG-specific enhancer or prevent in-vivo crosstalk.

## Reproduction

Run `.venv/bin/python scripts/prepare_cn4461_pilot_and_barcodes.py` from the
repository root. It freezes the unchanged selected pilot exports, regenerates
the codebook with a fixed seed, exhaustively checks distances and simulates all
30 cassette replacements. `barcode_replacement_design.json`, `codebook_QC.json`
and `engineering_simulation_QC.json` contain sequence-level results. No wet-lab
measurements are represented as completed.

Original architecture: [CN4461 reference](https://www.addgene.org/236600/sequences/)
and [Hunker et al., Neuron 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12237622/).
Enzyme cut definitions: [NEB](https://www.neb.com/en/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities).
