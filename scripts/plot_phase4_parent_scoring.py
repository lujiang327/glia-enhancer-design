#!/usr/bin/env python3
"""Plot compact parent-sequence prediction QC."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    table = pd.read_csv(args.predictions, sep="\t")
    observed = np.log2(table["muller_mean_cpm"].to_numpy() + 0.25)
    forward = table["chrombpnet_forward_logcount"].to_numpy()
    reverse = table["chrombpnet_reverse_complement_logcount"].to_numpy()
    mean = table["chrombpnet_mean_logcount"].to_numpy()
    delta = table["chrombpnet_orientation_abs_logcount_delta"].to_numpy()

    figure, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    scatter = axes[0].scatter(observed, mean, c=delta, cmap="viridis", s=13, alpha=0.75)
    axes[0].set_xlabel("Observed Muller accessibility, log2(CPM + 0.25)")
    axes[0].set_ylabel("Strand-mean ChromBPNet log-count")
    axes[0].set_title("Parent prediction vs observed")
    figure.colorbar(scatter, ax=axes[0], label="Absolute F/RC log-count delta")

    axes[1].scatter(forward, reverse, s=13, alpha=0.65, color="#4c78a8")
    low = min(forward.min(), reverse.min()); high = max(forward.max(), reverse.max())
    axes[1].plot([low, high], [low, high], linestyle="--", color="black", linewidth=1)
    axes[1].set_xlabel("Forward log-count")
    axes[1].set_ylabel("Reverse-complement log-count")
    axes[1].set_title("Orientation consistency")

    axes[2].hist(delta, bins=40, color="#f58518", edgecolor="white")
    axes[2].axvline(1, linestyle="--", color="black", linewidth=1, label="flag threshold")
    axes[2].set_xlabel("Absolute F/RC log-count delta")
    axes[2].set_ylabel("Candidates")
    axes[2].set_title("Orientation uncertainty")
    axes[2].legend(frameon=False)
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
