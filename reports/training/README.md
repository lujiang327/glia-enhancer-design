# Müller glia model-training QC

The fold 0, threshold 0.5, seed 42 Tn5 bias model is **accepted for the full
ChromBPNet pilot with a GC-composition caveat**. A fold 0 Müller-glia model has
now been trained with a 1.0 negative-sampling ratio and passed initial held-out
numerical and motif-plausibility screens. It remains a provisional candidate;
differential accessibility remains blocked pending seed reproducibility and
donor-sensitivity checks.

Held-out numerical QC on `chr1`, `chr3`, and `chr6` passed: nonpeak count
Pearson/Spearman correlations were 0.625/0.662, and peak count correlation was
positive (Pearson 0.359), avoiding the problematic anticorrelation failure mode.

TF-MoDISco used 30,000 deterministic held-out peak regions. All 15 profile
patterns had a Tn5 motif as their top database match; the first two contain
57.2% of profile seqlets and the first five contain 84.5%. The dominant learned
profile matrices visually match Tn5 patterns. Count attributions are more
composition-driven: 72.0% of positive count seqlets fall in patterns with an
enzyme-bias top match, while most remaining matrices are broad GC-rich,
AT-rich, or low-complexity effects.

The main caveat is a low-q SP2 match for one negative count pattern (778
seqlets, 1.97% of all count seqlets) and two small low-q ZN770 matches (165
seqlets together, 0.42%). Their learned matrices are diffuse composition
effects rather than compact versions of the matched TF logos. They do not block
the pilot, but the full no-bias model must be checked for residual Tn5 signal
and excessive GC-rich motifs. A systematic GC signal in its leading motifs
would trigger bias-threshold reconsideration.

See `bias_fold0_seed42_training.json`,
`bias_fold0_seed42_numerical_qc.json`, and
`bias_fold0_seed42_motif_qc.json` for machine-readable evidence.

Full-model preparation retained 256,225 peaks and 256,225 matched nonpeaks.
Only 23 of 256,248 prepared peaks were removed, the count-loss weight is 25.4,
and the saved parameters retain 500-bp jitter with the planned 512-filter,
eight-dilation architecture. See `chrombpnet_fold0_preparation.json`.

The balanced candidate trained for 15 epochs and restored epoch 10 (one-based),
with validation loss improving from 422.26 to 379.17. On held-out chromosomes
1, 3, and 6, peak count Pearson correlation was 0.668 and combined peak/nonpeak
correlation was 0.776. Its nonpeak count MSE was 1.659, a 49.4% improvement over
the initial 0.1-negative-ratio model, while peak performance declined modestly.
The full model still adds signal at some nonpeaks, so calibration remains a
caveat.

TF-MoDISco recovered NFI, SOX, LHX2/homeobox, AP-1, RFX, RORA, CRX/OTX,
CTCF, SP/KLF, and NFY families. No motif pattern had Tn5 as its top match.
Müller-glia RNA supports the NFI, SOX9, and LHX2 programs. CRX/OTX and generic
homeobox patterns require caution because they may reflect shared retinal
grammar or donor-specific contamination. See
`chrombpnet_fold0_seed42_neg1_motif_qc.json` and its companion Markdown report.

![Selected learned bias motifs and database matches](bias_motif_qc_selected.png)

Representative full-model contribution matrices are in
`chrombpnet_neg1_motif_logos/`.
