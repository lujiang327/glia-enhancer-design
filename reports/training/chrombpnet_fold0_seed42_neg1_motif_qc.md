# Müller-glia ChromBPNet motif QC: fold 0, seed 42, negative ratio 1

**Result: PASS with caveats.** The balanced model learned biologically plausible
retinal and glial sequence features, and its no-bias checkpoint shows no
dominant residual Tn5 motif. It is the provisional primary candidate for
reproducibility testing. This result does not yet open the differential
accessibility phase.

Attributions were calculated for 30,000 held-out peaks on chromosomes 1, 3,
and 6. Profile TF-MoDISco found 50 patterns from 28,182 seqlets; count
TF-MoDISco found 38 patterns from 22,191 seqlets. The leading contribution
matrices include NFI, SOX, LHX2/homeobox, AP-1, RFX, RORA/nuclear receptor,
CRX/OTX, CTCF, SP/KLF, and NFY families. Their compact matrices and repeated
recovery across profile and count tasks support learned biological sequence
grammar.

RNA provides strong Müller-glia support for `NFIA`, `NFIB`, `NFIX`, `SOX9`,
and `LHX2`; `RFX3`, `FOS`, `JUNB`, `RORA`, and `PAX6` are also detected. The
strong CRX/OTX and generic homeobox patterns need caution. `CRX` and `OTX2`
are less frequent in the published Müller-glia barcodes than in other retinal
cells and are elevated in donor LGS3. Matches named NKX6-1 or NOTO indicate
homeobox-family similarity because those genes are not detected in the
Müller-glia RNA data.

No pattern has Tn5 as its top match. A small negative profile pattern matched
DNASE_2 (41 seqlets; 0.145% of profile seqlets), and no count pattern matched
DNase or Tn5. This does not suggest residual enzyme bias dominates the learned
attributions.

The model still requires independent-seed consistency and donor-sensitivity
checks. Until those pass, its status remains candidate QC in progress and
differential accessibility remains blocked.

Representative contribution matrices are in
`chrombpnet_neg1_motif_logos/`; the motif-family RNA table is
`chrombpnet_neg1_motif_rna_support.tsv`.
