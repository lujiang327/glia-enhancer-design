#!/usr/bin/env python3
"""Nominate and directly score dual-orientation single-base perturbations.

The pilot and full scan share this implementation; scope comes from config.
"""

import argparse
import hashlib
import json
from pathlib import Path


BASES = "ACGT"


def reverse_complement_one_hot(values):
    return values[:, ::-1, ::-1]


def select_spaced(records, count, minimum_spacing):
    """Greedily select high-scoring edits at separated parent positions."""
    selected = []
    for record in sorted(records, key=lambda value: (-value["nomination_score"], value["parent_position_0based"])):
        if all(abs(record["parent_position_0based"] - prior["parent_position_0based"]) >= minimum_spacing
               for prior in selected):
            selected.append(record)
            if len(selected) == count:
                break
    return selected


def nominate_edits(sequence_one_hot, count_hypothetical, profile_hypothetical,
                    parent_start, parent_end, gain_count, loss_count, minimum_spacing):
    """Nominate count-gain substitutions and matched loss controls."""
    gain, loss = [], []
    for absolute_position in range(parent_start, parent_end):
        parent_position = absolute_position - parent_start
        ref_index = int(sequence_one_hot[absolute_position].argmax())
        if sequence_one_hot[absolute_position].sum() != 1:
            continue
        alternatives = [index for index in range(4) if index != ref_index]
        count_current = float(count_hypothetical[absolute_position, ref_index])
        profile_current = float(profile_hypothetical[absolute_position, ref_index])
        best = max(alternatives, key=lambda index: count_hypothetical[absolute_position, index])
        worst = min(alternatives, key=lambda index: count_hypothetical[absolute_position, index])
        gain_delta = float(count_hypothetical[absolute_position, best] - count_current)
        loss_delta = float(count_hypothetical[absolute_position, worst] - count_current)
        if gain_delta > 0:
            gain.append({
                "design_class": "gain", "parent_position_0based": parent_position,
                "ref": BASES[ref_index], "alt": BASES[best], "nomination_score": gain_delta,
                "nomination_count_hypothetical_delta": gain_delta,
                "nomination_profile_hypothetical_delta": float(
                    profile_hypothetical[absolute_position, best] - profile_current),
            })
        if loss_delta < 0:
            loss.append({
                "design_class": "loss_control", "parent_position_0based": parent_position,
                "ref": BASES[ref_index], "alt": BASES[worst], "nomination_score": -loss_delta,
                "nomination_count_hypothetical_delta": loss_delta,
                "nomination_profile_hypothetical_delta": float(
                    profile_hypothetical[absolute_position, worst] - profile_current),
            })
    return (select_spaced(gain, gain_count, minimum_spacing)
            + select_spaced(loss, loss_count, minimum_spacing))


def identify_outputs(model, predictions):
    if not isinstance(predictions, (list, tuple)) or len(predictions) != 2:
        raise ValueError("Expected two ChromBPNet outputs")
    names = list(model.output_names)
    shapes = [value.shape for value in predictions]
    profile_index = next((i for i, shape in enumerate(shapes) if len(shape) == 2 and shape[1] > 1), None)
    count_index = next((i for i, shape in enumerate(shapes) if len(shape) == 2 and shape[1] == 1), None)
    if profile_index is None or count_index is None:
        raise ValueError("Could not identify profile and count outputs: {}".format(names))
    return predictions[profile_index], predictions[count_index]


def softmax(logits):
    import numpy as np
    centered = logits - np.max(logits, axis=1, keepdims=True)
    values = np.exp(centered)
    return values / values.sum(axis=1, keepdims=True)


def jsd(left, right):
    from scipy import stats
    midpoint = (left + right) / 2
    return float(stats.entropy(midpoint) - (stats.entropy(left) + stats.entropy(right)) / 2)


def sha256_text(value):
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def edit_screen_pass(minimum_effect, forward_jsd, reverse_jsd, min_effect, max_jsd):
    return (minimum_effect >= min_effect) & (forward_jsd <= max_jsd) & (reverse_jsd <= max_jsd)


def load_pilot_candidates(qc_path, parents_path, count):
    import pandas as pd
    qc = pd.read_csv(qc_path, sep="\t")
    parents = pd.read_csv(parents_path, sep="\t")
    keys = ["peak_id", "candidate_rank", "chrom", "start", "end"]
    if qc["peak_id"].duplicated().any() or parents["peak_id"].duplicated().any():
        raise ValueError("Duplicate parent identifiers")
    candidates = qc.merge(
        parents[keys + ["parent_sequence", "parent_sequence_sha256"]],
        on=keys, how="left", validate="one_to_one",
    )
    if candidates["parent_sequence"].isna().any():
        raise ValueError("QC rows do not match parent identifiers and coordinates")
    for row in candidates.itertuples(index=False):
        if sha256_text(row.parent_sequence) != row.parent_sequence_sha256:
            raise ValueError("Parent sequence checksum mismatch for {}".format(row.peak_id))
    pilot = candidates[candidates["attribution_qc_pass"]].sort_values("candidate_rank").head(count).copy()
    if len(pilot) != count:
        raise ValueError("Insufficient attribution-passing parents for pilot")
    return pilot


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-qc", required=True, type=Path)
    parser.add_argument("--parents", required=True, type=Path)
    parser.add_argument("--attributions", required=True, type=Path)
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    import h5py
    import numpy as np
    import pandas as pd
    import tensorflow as tf
    import chrombpnet.training.utils.losses as losses
    from tensorflow.keras.models import load_model
    from tensorflow.keras.utils import get_custom_objects

    config = json.loads(args.config.read_text())
    parent_count = int(config.get("parent_count", config.get("pilot_parents", 100)))
    pilot = load_pilot_candidates(args.candidate_qc, args.parents, parent_count)

    with h5py.File(args.attributions, "r") as handle:
        h5_ids = [value.decode() if isinstance(value, bytes) else str(value) for value in handle["peak_id"][:]]
        h5_index = {value: index for index, value in enumerate(h5_ids)}
        indices = np.asarray([h5_index[value] for value in pilot["peak_id"]], dtype=int)
        if len(h5_index) != len(h5_ids):
            raise ValueError("Duplicate attribution identifiers")
        sequences = handle["sequence_one_hot"][indices].astype(np.int8)
        count_hypothetical = handle["counts/consensus_hypothetical"][indices].astype(np.float32)
        profile_hypothetical = handle["profile/consensus_hypothetical"][indices].astype(np.float32)
    if not np.all(sequences.sum(axis=2) == 1):
        raise ValueError("Selected model inputs contain ambiguous bases")
    if not np.isfinite(count_hypothetical).all() or not np.isfinite(profile_hypothetical).all():
        raise ValueError("Non-finite hypothetical attribution")

    parent_start, parent_end = config["parent_slice_zero_based_half_open"]
    edit_rows = []
    model_inputs = []
    sequence_ids = []
    parent_input_index = {}
    nomination_rows = []
    for local_index, row in enumerate(pilot.itertuples(index=False)):
        parent_one_hot = sequences[local_index]
        parent_sequence = "".join(BASES[index] for index in parent_one_hot[parent_start:parent_end].argmax(axis=1))
        if parent_sequence != row.parent_sequence:
            raise ValueError("Attribution sequence mismatch for {}".format(row.peak_id))
        parent_input_index[row.peak_id] = len(model_inputs)
        model_inputs.append(parent_one_hot)
        sequence_ids.append("{}:parent".format(row.peak_id))
        edits = nominate_edits(
            parent_one_hot, count_hypothetical[local_index], profile_hypothetical[local_index],
            parent_start, parent_end, int(config["gain_substitutions_per_parent"]),
            int(config["loss_control_substitutions_per_parent"]),
            int(config["minimum_spacing_bp_within_class"]),
        )
        expected = int(config["gain_substitutions_per_parent"]) + int(config["loss_control_substitutions_per_parent"])
        if len(edits) != expected and not config.get("allow_incomplete_nominations", False):
            raise ValueError("Could not nominate {} edits for {}".format(expected, row.peak_id))
        nomination_rows.append({
            "candidate_rank": row.candidate_rank, "peak_id": row.peak_id,
            "gain_nominated": sum(edit["design_class"] == "gain" for edit in edits),
            "loss_control_nominated": sum(edit["design_class"] == "loss_control" for edit in edits),
            "nomination_complete": len(edits) == expected,
        })
        for edit_number, edit in enumerate(edits, 1):
            position = int(edit["parent_position_0based"])
            if parent_sequence[position] != edit["ref"]:
                raise ValueError("Reference allele mismatch")
            variant_sequence = parent_sequence[:position] + edit["alt"] + parent_sequence[position + 1:]
            variant_input = parent_one_hot.copy()
            variant_input[parent_start + position, :] = 0
            variant_input[parent_start + position, BASES.index(edit["alt"])] = 1
            variant_id = "{}:{}:{}{}>{}".format(
                row.peak_id, edit["design_class"], position, edit["ref"], edit["alt"]
            )
            variant_input_index = len(model_inputs)
            model_inputs.append(variant_input)
            sequence_ids.append(variant_id)
            edit_rows.append({
                "candidate_rank": row.candidate_rank, "peak_id": row.peak_id,
                "chrom": row.chrom, "start": row.start, "end": row.end,
                "genomic_class": row.genomic_class, "ranking_score": row.ranking_score,
                "muller_mean_cpm": row.muller_mean_cpm,
                "observed_specificity_log2_ratio": row.observed_specificity_log2_ratio,
                "parent_count_orientation_flag": row.parent_count_orientation_flag,
                "variant_off_target_specificity_prediction_available": False,
                "parent_sequence": parent_sequence, "parent_sequence_sha256": sha256_text(parent_sequence),
                "variant_id": variant_id, "variant_sequence": variant_sequence,
                "variant_sequence_sha256": sha256_text(variant_sequence),
                "edit_number_within_parent": edit_number,
                "genomic_position_0based": int(row.start) + position,
                "parent_input_index": parent_input_index[row.peak_id],
                "variant_input_index": variant_input_index,
                **edit,
            })

    inputs = np.asarray(model_inputs, dtype=np.int8)
    get_custom_objects().update({"multinomial_nll": losses.multinomial_nll, "tf": tf})
    model = load_model(str(args.model), compile=False)
    if inputs.shape[1] != int(model.input_shape[1]):
        raise ValueError("Model input length mismatch")
    batch_size = int(config["batch_size"])
    forward_raw = model.predict(inputs, batch_size=batch_size, verbose=1)
    reverse_raw = model.predict(reverse_complement_one_hot(inputs), batch_size=batch_size, verbose=1)
    forward_logits, forward_counts = identify_outputs(model, forward_raw)
    reverse_logits, reverse_counts = identify_outputs(model, reverse_raw)
    forward_profiles = softmax(np.asarray(forward_logits, dtype=np.float64))
    reverse_profiles = softmax(np.asarray(reverse_logits, dtype=np.float64))[:, ::-1]
    forward_counts = np.asarray(forward_counts, dtype=float).reshape(-1)
    reverse_counts = np.asarray(reverse_counts, dtype=float).reshape(-1)
    if not all(np.all(np.isfinite(value)) for value in (
        forward_profiles, reverse_profiles, forward_counts, reverse_counts
    )):
        raise ValueError("Non-finite model predictions")
    for row in pilot.itertuples(index=False):
        i = parent_input_index[row.peak_id]
        if not np.isclose((forward_counts[i] + reverse_counts[i]) / 2,
                          row.chrombpnet_mean_logcount, atol=1e-5, rtol=0):
            raise ValueError("Rescored parent differs from frozen baseline for {}".format(row.peak_id))

    if not edit_rows:
        raise ValueError("No substitutions nominated; preserve failed run for review")

    results = pd.DataFrame(edit_rows)
    metrics = []
    for row in results.itertuples(index=False):
        parent_index, variant_index = int(row.parent_input_index), int(row.variant_input_index)
        delta_forward = forward_counts[variant_index] - forward_counts[parent_index]
        delta_reverse = reverse_counts[variant_index] - reverse_counts[parent_index]
        intended = 1 if row.design_class == "gain" else -1
        metrics.append({
            "parent_forward_logcount": forward_counts[parent_index],
            "parent_reverse_complement_logcount": reverse_counts[parent_index],
            "parent_mean_logcount": (forward_counts[parent_index] + reverse_counts[parent_index]) / 2,
            "variant_forward_logcount": forward_counts[variant_index],
            "variant_reverse_complement_logcount": reverse_counts[variant_index],
            "variant_mean_logcount": (forward_counts[variant_index] + reverse_counts[variant_index]) / 2,
            "delta_forward_logcount": delta_forward,
            "delta_reverse_complement_logcount": delta_reverse,
            "delta_mean_logcount": (delta_forward + delta_reverse) / 2,
            "robust_effect_lower_bound": min(intended * delta_forward, intended * delta_reverse),
            "intended_effect_in_both_orientations": bool(
                intended * delta_forward > 0 and intended * delta_reverse > 0),
            "orientation_effect_abs_difference": abs(delta_forward - delta_reverse),
            "forward_profile_jsd_vs_parent": jsd(
                forward_profiles[parent_index], forward_profiles[variant_index]),
            "reverse_profile_jsd_vs_parent": jsd(
                reverse_profiles[parent_index], reverse_profiles[variant_index]),
        })
    results = pd.concat([results.reset_index(drop=True), pd.DataFrame(metrics)], axis=1)
    results["mean_profile_jsd_vs_parent"] = (
        results["forward_profile_jsd_vs_parent"] + results["reverse_profile_jsd_vs_parent"]
    ) / 2
    if "edit_screen" in config:
        screen = config["edit_screen"]
        results["edit_screen_pass"] = edit_screen_pass(
            results["robust_effect_lower_bound"], results["forward_profile_jsd_vs_parent"],
            results["reverse_profile_jsd_vs_parent"],
            float(screen["minimum_intended_logcount_effect_each_orientation"]),
            float(screen["maximum_profile_jsd_vs_parent_each_orientation"]),
        )
    results = results.sort_values(
        ["design_class", "robust_effect_lower_bound", "candidate_rank"], ascending=[True, False, True]
    )

    args.output_dir.mkdir(parents=True, exist_ok=False)
    prefix = config.get("output_prefix", "perturbation_pilot")
    table_path = args.output_dir / (prefix + ".tsv.gz")
    results.to_csv(table_path, sep="\t", index=False, compression="gzip")
    profile_path = args.output_dir / (prefix + "_profiles.h5")
    string_type = h5py.string_dtype(encoding="utf-8")
    with h5py.File(profile_path, "w") as handle:
        handle.create_dataset("sequence_id", data=np.asarray(sequence_ids, dtype=object), dtype=string_type)
        handle.create_dataset("forward_profile_probability", data=forward_profiles.astype(np.float32), compression="gzip")
        handle.create_dataset("reverse_complement_profile_probability", data=reverse_profiles.astype(np.float32), compression="gzip")

    by_class = {}
    for design_class, group in results.groupby("design_class"):
        by_class[design_class] = {
            "variants": int(len(group)),
            "intended_effect_both_orientations": int(group["intended_effect_in_both_orientations"].sum()),
            "intended_effect_fraction": float(group["intended_effect_in_both_orientations"].mean()),
            "median_delta_mean_logcount": float(group["delta_mean_logcount"].median()),
            "median_robust_effect_lower_bound": float(group["robust_effect_lower_bound"].median()),
            "median_profile_jsd_vs_parent": float(group["mean_profile_jsd_vs_parent"].median()),
        }
        if "edit_screen_pass" in group:
            passing = group[group["edit_screen_pass"]]
            by_class[design_class]["screen_pass_variants"] = int(len(passing))
            by_class[design_class]["screen_pass_parents"] = int(passing["peak_id"].nunique())
    coverage = pd.DataFrame(nomination_rows)
    if "edit_screen_pass" in results:
        for design_class in ("gain", "loss_control"):
            passing = results[(results["design_class"] == design_class) & results["edit_screen_pass"]]
            coverage[design_class + "_screen_pass"] = coverage["peak_id"].map(passing.groupby("peak_id").size()).fillna(0).astype(int)
    coverage_path = args.output_dir / (prefix + "_parent_coverage.tsv.gz")
    coverage.to_csv(coverage_path, sep="\t", index=False, compression="gzip")
    summary = {
        "status": config.get("result_status", "PERTURBATION_PILOT_READY_FOR_REVIEW"),
        "parents": int(len(pilot)), "variants": int(len(results)),
        "model_sequences_scored_per_orientation": int(len(inputs)),
        "by_design_class": by_class,
        "parents_with_incomplete_nominations": int((~coverage["nomination_complete"]).sum()),
        "outputs": {"table": str(table_path), "profiles": str(profile_path), "coverage": str(coverage_path)},
        "next_gate": config["next_gate"],
    }
    (args.output_dir / (prefix + "_summary.json")).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
