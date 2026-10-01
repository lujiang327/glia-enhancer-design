#!/usr/bin/env python3
"""Extract natural candidate sequences and score them with frozen ChromBPNet."""

import argparse
import hashlib
import json
import math
from pathlib import Path


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def sha256_text(value):
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def quantiles(values, probabilities):
    import numpy as np
    return {str(probability): float(np.quantile(values, probability)) for probability in probabilities}


def identify_outputs(model, predictions):
    if not isinstance(predictions, (list, tuple)) or len(predictions) != 2:
        raise ValueError("Expected two ChromBPNet outputs")
    names = list(model.output_names)
    by_name = dict(zip(names, predictions))
    profile_name = next((name for name in names if "profile" in name), None)
    count_name = next((name for name in names if "count" in name), None)
    if profile_name is None or count_name is None or profile_name == count_name:
        shapes = [value.shape for value in predictions]
        profile_index = next((index for index, shape in enumerate(shapes) if len(shape) == 2 and shape[1] > 1), None)
        count_index = next((index for index, shape in enumerate(shapes) if len(shape) == 2 and shape[1] == 1), None)
        if profile_index is None or count_index is None:
            raise ValueError("Could not identify profile and count outputs")
        return predictions[profile_index], predictions[count_index], names
    return by_name[profile_name], by_name[count_name], names


def softmax(logits):
    import numpy as np
    centered = logits - np.max(logits, axis=1, keepdims=True)
    exp = np.exp(centered)
    return exp / np.sum(exp, axis=1, keepdims=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    import h5py
    import numpy as np
    import pandas as pd
    import pyfaidx
    from scipy import stats
    import tensorflow as tf
    import chrombpnet.training.utils.losses as losses
    from chrombpnet.training.utils.one_hot import dna_to_one_hot
    from tensorflow.keras.models import load_model
    from tensorflow.keras.utils import get_custom_objects

    config = json.loads(args.config.read_text())
    candidates = pd.read_csv(args.candidates, sep="\t")
    expected = int(config["parent_sequence_bp"])
    if len(candidates) != 1000 or candidates["peak_id"].nunique() != len(candidates):
        raise ValueError("Expected 1,000 unique candidate rows")
    if not np.array_equal(candidates["candidate_rank"].to_numpy(), np.arange(1, len(candidates) + 1)):
        raise ValueError("Candidate ranks must be exactly 1..1000")
    if np.any(candidates["end"].to_numpy() - candidates["start"].to_numpy() != expected):
        raise ValueError("Candidate parent intervals must be 500 bp")

    get_custom_objects().update({"multinomial_nll": losses.multinomial_nll, "tf": tf})
    model = load_model(str(args.model), compile=False)
    input_length = int(model.input_shape[1])
    if input_length != int(config["model_input_bp"]):
        raise ValueError("Unexpected model input length: {}".format(input_length))

    genome = pyfaidx.Fasta(str(args.genome))
    parent_sequences = []
    model_sequences = []
    input_starts = []
    input_ends = []
    invalid_parent_bases = []
    invalid_input_bases = []
    for row in candidates.itertuples(index=False):
        center = int(row.start) + expected // 2
        input_start = center - input_length // 2
        input_end = center + input_length // 2
        parent = str(genome[row.chrom][int(row.start):int(row.end)]).upper()
        model_sequence = str(genome[row.chrom][input_start:input_end]).upper()
        if len(parent) != expected or len(model_sequence) != input_length:
            raise ValueError("Reference extraction failed for {}".format(row.peak_id))
        parent_sequences.append(parent)
        model_sequences.append(model_sequence)
        input_starts.append(input_start)
        input_ends.append(input_end)
        invalid_parent_bases.append(sum(base not in "ACGT" for base in parent))
        invalid_input_bases.append(sum(base not in "ACGT" for base in model_sequence))
    genome.close()

    forward = dna_to_one_hot(model_sequences)
    reverse = forward[:, ::-1, ::-1]
    batch_size = int(config["batch_size"])
    forward_predictions = model.predict(forward, batch_size=batch_size, verbose=1)
    reverse_predictions = model.predict(reverse, batch_size=batch_size, verbose=1)
    forward_logits, forward_logcounts, output_names = identify_outputs(model, forward_predictions)
    reverse_logits, reverse_logcounts, reverse_names = identify_outputs(model, reverse_predictions)
    if output_names != reverse_names:
        raise ValueError("Model output names changed between orientations")

    forward_profiles = softmax(np.asarray(forward_logits, dtype=np.float64))
    reverse_profiles = softmax(np.asarray(reverse_logits, dtype=np.float64))[:, ::-1]
    mean_profiles = (forward_profiles + reverse_profiles) / 2
    forward_logcounts = np.asarray(forward_logcounts, dtype=np.float64).reshape(-1)
    reverse_logcounts = np.asarray(reverse_logcounts, dtype=np.float64).reshape(-1)
    mean_logcounts = (forward_logcounts + reverse_logcounts) / 2
    if not all(np.all(np.isfinite(value)) for value in (forward_profiles, reverse_profiles, mean_profiles,
                                                        forward_logcounts, reverse_logcounts, mean_logcounts)):
        raise ValueError("Non-finite model prediction")

    observed = np.log2(candidates["muller_mean_cpm"].to_numpy(dtype=float) + 0.25)
    pearson = float(stats.pearsonr(mean_logcounts, observed)[0])
    spearman = float(stats.spearmanr(mean_logcounts, observed)[0])
    orientation_pearson = float(stats.pearsonr(forward_logcounts, reverse_logcounts)[0])
    orientation_abs_delta = np.abs(forward_logcounts - reverse_logcounts)
    profile_orientation_jsd = []
    for forward_profile, reverse_profile in zip(forward_profiles, reverse_profiles):
        profile_orientation_jsd.append(float(stats.entropy((forward_profile + reverse_profile) / 2)
                                             - (stats.entropy(forward_profile) + stats.entropy(reverse_profile)) / 2))
    profile_orientation_jsd = np.asarray(profile_orientation_jsd)

    output = candidates.copy()
    output["parent_sequence"] = parent_sequences
    output["parent_sequence_sha256"] = [sha256_text(sequence) for sequence in parent_sequences]
    output["model_input_start"] = input_starts
    output["model_input_end"] = input_ends
    output["parent_invalid_bases"] = invalid_parent_bases
    output["model_input_invalid_bases"] = invalid_input_bases
    output["chrombpnet_forward_logcount"] = forward_logcounts
    output["chrombpnet_reverse_complement_logcount"] = reverse_logcounts
    output["chrombpnet_mean_logcount"] = mean_logcounts
    output["chrombpnet_orientation_abs_logcount_delta"] = orientation_abs_delta
    output["chrombpnet_orientation_profile_jsd"] = profile_orientation_jsd

    args.output_dir.mkdir(parents=True, exist_ok=False)
    table_path = args.output_dir / "parent_predictions.tsv.gz"
    output.to_csv(table_path, sep="\t", index=False, compression="gzip")
    profile_path = args.output_dir / "parent_profiles.h5"
    string_type = h5py.string_dtype(encoding="utf-8")
    with h5py.File(profile_path, "w") as handle:
        handle.create_dataset("peak_id", data=candidates["peak_id"].astype(str).to_numpy(), dtype=string_type)
        handle.create_dataset("candidate_rank", data=candidates["candidate_rank"].to_numpy(dtype=np.int32))
        handle.create_dataset("forward_profile_probability", data=forward_profiles.astype(np.float32), compression="gzip")
        handle.create_dataset("reverse_complement_profile_probability", data=reverse_profiles.astype(np.float32), compression="gzip")
        handle.create_dataset("mean_profile_probability", data=mean_profiles.astype(np.float32), compression="gzip")

    probabilities = (0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)
    summary = {
        "status": "PARENT_PREDICTIONS_READY_FOR_QC",
        "candidates": len(candidates),
        "model_input_length": input_length,
        "model_output_length": int(mean_profiles.shape[1]),
        "model_output_names": output_names,
        "parent_sequence_length": expected,
        "parents_with_invalid_bases": int(np.sum(np.asarray(invalid_parent_bases) > 0)),
        "inputs_with_invalid_bases": int(np.sum(np.asarray(invalid_input_bases) > 0)),
        "total_parent_invalid_bases": int(np.sum(invalid_parent_bases)),
        "total_input_invalid_bases": int(np.sum(invalid_input_bases)),
        "prediction_vs_observed_muller_accessibility": {
            "pearsonr": pearson,
            "spearmanr": spearman,
            "scope": "Association screen on selected candidates, not held-out model evaluation."
        },
        "reverse_complement_consistency": {
            "logcount_pearsonr": orientation_pearson,
            "absolute_logcount_delta_quantiles": quantiles(orientation_abs_delta, probabilities),
            "profile_jsd_quantiles": quantiles(profile_orientation_jsd, probabilities)
        },
        "mean_logcount_quantiles": quantiles(mean_logcounts, probabilities),
        "outputs": {"table": str(table_path), "profiles": str(profile_path)},
        "next_gate": config["next_gate"]
    }
    summary_path = args.output_dir / "parent_scoring_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
