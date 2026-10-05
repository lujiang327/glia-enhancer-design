# Motif-window attribution assessment

This CPU job connects the 59,828 edit-overlapping motif comparisons to the existing count/profile DeepSHAP attribution for all 1,000 natural parents. It preserves the 2,884-row library and adds edited-base attribution to its 1,884 variants. No model training, new differential testing or variant attribution is performed.

The job verifies SHA256 for the candidate library, motif comparison table and saved HDF5 attribution. It checks parent identifiers, ranks, sequences and the 500bp slice [807,1307) in the 2,114bp model input. Saved reverse-complement scores are already aligned to forward positions and A/C/G/T channels; they must not be reversed again.

Natural-base projection selects the observed base's hypothetical score at each position, separately for count and profile targets. Each motif window receives the sum of projected attribution for forward, aligned reverse and consensus maps, its share of the parent's absolute consensus attribution, and its overlap with the top 10% of parent bases by absolute attribution. Top-base ties use stable position order. A zero attribution denominator yields missing rather than invented fractional evidence.

Every edited base receives its parent projected score and absolute-magnitude percentile, plus alternate-minus-reference hypothetical attribution in both orientations and consensus. Hypothetical count-delta agreement with the directly predicted edit direction is recorded descriptively. The profile head describes shape; its sign is not an accessibility gain/loss criterion. Parent attribution and saved float16 hypothetical deltas are not recomputed variant SHAP or exact direct-model edit effects.

A gained motif window is evaluated on the natural parent sequence in that same interval. This reports pre-edit context, not attribution to the gained motif in the variant. Motif windows overlap and redundant motif models share bases: do not sum their attribution as independent regulatory contributions. Positive/negative window sums and top-base overlap are descriptive evidence, not a significance test or causal TF identification. No new hard cutoff or automatic sequence exclusion is introduced.

Outputs:

- `edit_motif_parent_attribution.tsv.gz`: the original motif comparisons with count/profile window attribution metrics.
- `ranked_library_with_parent_attribution.tsv.gz`: all original library evidence plus edited-base attribution metrics.
- `motif_attribution_summary.json`, runtime provenance and file checksums.

The account and scratch root remain manual environment settings. Submit from the Great Lakes repository root after pulling:

```bash
sbatch --account="$SLURM_ACCOUNT" hpc/slurm/annotate_phase5_motif_attribution.sbatch
```

`PHASE3_SCRATCH_ROOT` must be exported. The job requests one CPU, 8GB RAM, two hours maximum, no GPU; outputs go to `phase5_motif_attribution` under scratch. Existing output directories are never overwritten.

Afterward, review representative gain/control pairs for agreement between direct edit effects, attribution and motif changes, then propose experimental subsets with their natural parents. Off-target variant specificity and final reporter/oligo design remain pending.
