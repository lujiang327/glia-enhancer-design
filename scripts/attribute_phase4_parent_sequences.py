#!/usr/bin/env python3
"""Compute aligned forward/RC ChromBPNet attribution maps for parent sequences."""

import argparse
import json
from pathlib import Path


def reverse_complement_one_hot(values):
    return values[:, ::-1, ::-1]


def align_reverse_complement_scores(values):
    return values[:, ::-1, ::-1]


def row_agreement(forward, reverse, start, end):
    """Agreement of projected, basewise scores in the parent slice."""
    import numpy as np
    a = np.sum(forward[:, start:end, :], axis=2)
    b = np.sum(reverse[:, start:end, :], axis=2)
    pearson = np.full(a.shape[0], np.nan, dtype=float)
    cosine = np.full(a.shape[0], np.nan, dtype=float)
    sign_concordance = np.full(a.shape[0], np.nan, dtype=float)
    for index, (left, right) in enumerate(zip(a, b)):
        if np.std(left) > 0 and np.std(right) > 0:
            pearson[index] = np.corrcoef(left, right)[0, 1]
        denominator = np.linalg.norm(left) * np.linalg.norm(right)
        if denominator > 0:
            cosine[index] = np.dot(left, right) / denominator
        magnitude = np.maximum(np.abs(left), np.abs(right))
        threshold = np.quantile(magnitude, 0.9)
        active = magnitude >= threshold
        if np.any(active):
            sign_concordance[index] = np.mean(np.sign(left[active]) == np.sign(right[active]))
    return pearson, cosine, sign_concordance


def finite_quantiles(values):
    import numpy as np
    probabilities = (0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)
    finite = np.asarray(values)[np.isfinite(values)]
    return {
        "finite": int(len(finite)),
        "quantiles": {str(p): float(np.quantile(finite, p)) for p in probabilities} if len(finite) else {},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    import deepdish as dd  # noqa: F401 - validates the pinned attribution runtime
    import h5py
    import numpy as np
    import pandas as pd
    import pyfaidx
    import shap
    import tensorflow as tf
    import chrombpnet.evaluation.interpret.shap_utils as shap_utils
    import chrombpnet.training.utils.losses as losses
    from chrombpnet.training.utils.one_hot import dna_to_one_hot
    from tensorflow.keras.models import load_model
    from tensorflow.keras.utils import get_custom_objects

    tf.compat.v1.disable_eager_execution()
    config = json.loads(args.config.read_text())
    candidates = pd.read_csv(args.candidates, sep="\t")
    if len(candidates) != 1000 or candidates["peak_id"].nunique() != 1000:
        raise ValueError("Expected 1,000 unique parent candidates")
    input_length = int(config["model_input_bp"])
    parent_start, parent_end = config["parent_slice_zero_based_half_open"]
    if parent_end - parent_start != int(config["parent_sequence_bp"]):
        raise ValueError("Parent slice does not match parent sequence length")

    genome = pyfaidx.Fasta(str(args.genome))
    sequences = []
    for row in candidates.itertuples(index=False):
        center = int(row.start) + int(config["parent_sequence_bp"]) // 2
        sequence = str(genome[row.chrom][center - input_length // 2:center + input_length // 2]).upper()
        if len(sequence) != input_length:
            raise ValueError("Reference extraction failed for {}".format(row.peak_id))
        observed_parent = sequence[parent_start:parent_end]
        if observed_parent != row.parent_sequence:
            raise ValueError("Parent sequence mismatch for {}".format(row.peak_id))
        sequences.append(sequence)
    genome.close()
    forward_sequences = dna_to_one_hot(sequences)
    reverse_sequences = reverse_complement_one_hot(forward_sequences)

    get_custom_objects().update({"multinomial_nll": losses.multinomial_nll, "tf": tf})
    model = load_model(str(args.model), compile=False)
    if int(model.input_shape[1]) != input_length:
        raise ValueError("Model input length mismatch")

    count_target = tf.reduce_sum(model.outputs[1], axis=-1)
    profile_target = shap_utils.get_weightedsum_meannormed_logits(model)

    def explain_pair(target, label):
        print("Building {} explainer".format(label), flush=True)
        explainer = shap.explainers.deep.TFDeepExplainer(
            (model.input, target),
            shap_utils.shuffle_several_times,
            combine_mult_and_diffref=shap_utils.combine_mult_and_diffref,
        )
        results = []
        for orientation, values in (("forward", forward_sequences), ("reverse complement", reverse_sequences)):
            print("Scoring {} {}".format(label, orientation), flush=True)
            result = np.asarray(explainer.shap_values(values, progress_message=100), dtype=np.float32)
            if result.shape != values.shape or not np.all(np.isfinite(result)):
                raise ValueError("Invalid attribution output for {} {}".format(label, orientation))
            results.append(result)
        return results

    count_forward, count_reverse = explain_pair(count_target, "count")
    count_reverse_aligned = align_reverse_complement_scores(count_reverse)
    profile_forward, profile_reverse = explain_pair(profile_target, "profile")
    profile_reverse_aligned = align_reverse_complement_scores(profile_reverse)
    count_consensus = (count_forward + count_reverse_aligned) / 2
    profile_consensus = (profile_forward + profile_reverse_aligned) / 2
    count_projected_forward = count_forward * forward_sequences
    count_projected_reverse = count_reverse_aligned * forward_sequences
    profile_projected_forward = profile_forward * forward_sequences
    profile_projected_reverse = profile_reverse_aligned * forward_sequences

    count_pearson, count_cosine, count_sign = row_agreement(
        count_projected_forward, count_projected_reverse, parent_start, parent_end
    )
    profile_pearson, profile_cosine, profile_sign = row_agreement(
        profile_projected_forward, profile_projected_reverse, parent_start, parent_end
    )

    output = candidates[[
        "candidate_rank", "peak_id", "chrom", "start", "end", "genomic_class",
        "ranking_score", "muller_mean_cpm", "observed_specificity_log2_ratio",
        "chrombpnet_mean_logcount", "chrombpnet_orientation_abs_logcount_delta",
        "model_input_invalid_bases",
    ]].copy()
    output["count_attribution_parent_pearson"] = count_pearson
    output["count_attribution_parent_cosine"] = count_cosine
    output["count_attribution_parent_top10pct_sign_concordance"] = count_sign
    output["profile_attribution_parent_pearson"] = profile_pearson
    output["profile_attribution_parent_cosine"] = profile_cosine
    output["profile_attribution_parent_top10pct_sign_concordance"] = profile_sign
    output["parent_count_orientation_flag"] = (
        output["chrombpnet_orientation_abs_logcount_delta"]
        > float(config["parent_count_orientation_flag_threshold"])
    )

    args.output_dir.mkdir(parents=True, exist_ok=False)
    metrics_path = args.output_dir / "attribution_agreement.tsv.gz"
    output.to_csv(metrics_path, sep="\t", index=False, compression="gzip")
    h5_path = args.output_dir / "dual_orientation_attributions.h5"
    string_type = h5py.string_dtype(encoding="utf-8")
    with h5py.File(h5_path, "w") as handle:
        handle.attrs["model_input_bp"] = input_length
        handle.attrs["parent_start"] = parent_start
        handle.attrs["parent_end"] = parent_end
        handle.create_dataset("peak_id", data=candidates["peak_id"].astype(str).to_numpy(), dtype=string_type)
        handle.create_dataset("candidate_rank", data=candidates["candidate_rank"].to_numpy(dtype=np.int32))
        handle.create_dataset("sequence_one_hot", data=forward_sequences.astype(np.int8), compression="gzip")
        for head, forward, reverse, consensus in (
            ("counts", count_forward, count_reverse_aligned, count_consensus),
            ("profile", profile_forward, profile_reverse_aligned, profile_consensus),
        ):
            group = handle.create_group(head)
            group.create_dataset("forward_hypothetical", data=forward.astype(np.float16), compression="gzip")
            group.create_dataset("reverse_complement_aligned_hypothetical", data=reverse.astype(np.float16), compression="gzip")
            group.create_dataset("consensus_hypothetical", data=consensus.astype(np.float16), compression="gzip")

    summary = {
        "status": "DUAL_ORIENTATION_ATTRIBUTIONS_READY_FOR_QC",
        "candidates": len(candidates),
        "model_input_bp": input_length,
        "parent_slice": [parent_start, parent_end],
        "count_attribution_parent_agreement": {
            "pearson": finite_quantiles(count_pearson),
            "cosine": finite_quantiles(count_cosine),
            "top10pct_sign_concordance": finite_quantiles(count_sign),
        },
        "profile_attribution_parent_agreement": {
            "pearson": finite_quantiles(profile_pearson),
            "cosine": finite_quantiles(profile_cosine),
            "top10pct_sign_concordance": finite_quantiles(profile_sign),
        },
        "parent_count_orientation_flagged": int(output["parent_count_orientation_flag"].sum()),
        "outputs": {"attributions": str(h5_path), "agreement": str(metrics_path)},
        "next_gate": config["next_gate"],
    }
    summary_path = args.output_dir / "attribution_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
