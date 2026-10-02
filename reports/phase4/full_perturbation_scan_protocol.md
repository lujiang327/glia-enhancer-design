# Full Muller single-base perturbation scan

## Scope

Apply the successful pilot method to all 947 attribution-QC-passing parents.
The 53 failed parents remain in the original 1,000-candidate audit table.
The frozen no-bias model and consensus hypothetical attribution HDF5 are
unchanged. The first 100 parents are included again to check reproducibility
against the pilot results.

Nominate up to five gain substitutions and five loss controls per parent,
spacing nominated positions within each class by at least 5 bp. Each variant
contains one substitution within the original 500-bp parent and preserves the
surrounding model context. This is a selected substitution scan, not exhaustive
in-silico saturation mutagenesis.

## Screening and outputs

Directly score all parents and variants in both orientations. A gain passes
when both count changes are at least +0.1 log units. A loss control passes when
both are at most -0.1. For either class, profile JSD versus its parent must be
at most 0.05 in each orientation. These project screening thresholds were
selected during pilot review.

Preserve failed edits alongside passing edits. Report each parent's nomination
counts, whether nominations were complete, and passing gain/control counts.
Flag observed parent orientation sensitivity and explicitly mark variant
off-target prediction as unavailable. Verify rescored mean parent predictions
against the recorded frozen baseline within absolute tolerance 1e-5.

The scan can score at most 10,417 sequences per orientation: 947 parents and
9,470 variants. Outputs under the private scratch root's
`phase4_perturbation_full/` are:

- `perturbation_full.tsv.gz`: all nominated variants, parent sequences,
  alleles, predictions, effect and profile metrics, and screen results;
- `perturbation_full_parent_coverage.tsv.gz`: all 947 selected parents and
  their nomination and passing-edit counts;
- `perturbation_full_profiles.h5`: parent and variant profiles in both
  orientations, indexed by sequence ID;
- `perturbation_full_summary.json`, `runtime.txt`, and
  `output_checksums.sha256`: results and provenance.

The job requests one A40, 8 CPUs, 48 GB RAM, and a 2-hour limit. The completed
pilot took 57 seconds. Queue time and full-scan runtime remain uncertain.

## Review after completion

Verify checksums, finite predictions, one-base edits, parent hashes, and
coverage of all 947 parents. Compare the first 100 parents' effects with the
pilot; review passing-edit coverage, low-prediction parents, failed nominations,
profile disruptions, and effect disagreement between orientations. Do not
infer improved retinal specificity from a gain in the Muller model. Retain
Phase 3 observed parent specificity as supporting data with its original scope.

Repeat and mappability annotation remain required before library selection.
Further multi-base designs require separate direct rescoring because edit
effects may interact. Experimental activity and variant specificity remain
untested.
