# Five pilot candidate peak tracks

P1, P47, P140, P25 and P779 are shown in that order. These are the original
candidate ranks, not a new ranking.

Open `output/pdf/MG_five_pilot_candidate_peak_tracks.pdf` from the repository root.
The same pages are saved as PNG files in this folder.

1. **Pooled retinal tracks:** 2.5 kb around each parent, all13 cell types. Compare
   strength across columns and off-target signal down each column; all panels
   share one y-scale.
2. **Exact500-bp zoom:** focuses on the sequences to be cloned. This page has its
   own shared y-scale, explicitly stated; do not compare its heights directly to
   page1 without accounting for the scale.
3. **Donor glial tracks:** each donor's MG, astrocyte and microglia profiles. All
   panels share one donor-page scale. Compare within each donor and candidate.

Shading marks the cloned500 bp. Signal outside that interval is genomic context,
not included in the insert. Tracks retain the original global IPM/base scaling
and100-bp bins, without smoothing or individual autoscaling. Header CPM comes
from equal-donor TMM accessibility, a different normalization. Neither measure
establishes reporter activity, delivery or experimental cell-type specificity.

`five_candidate_metrics.tsv` records metrics; `five_candidate_track_means.tsv`
records pooled/donor means; `displayed_track_values.npz` preserves exact plotted
samples. The BED file provides the five hg38 intervals. Source hashes, settings
and cross-checks against the earlier audit are in `provenance_and_validation.json`.

Reproduce with `MPLCONFIGDIR=/tmp/glia-matplotlib .venv/bin/python scripts/plot_pilot_five_peak_tracks.py`.
