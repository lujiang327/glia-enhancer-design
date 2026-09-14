#!/usr/bin/env python3
"""Fail on non-finite Keras epoch metrics and summarize validation loss."""
import argparse
import csv
import json
import math
from pathlib import Path


def run(history_path, output_path):
    with open(history_path) as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError('Training history is empty')
    required = {'epoch', 'loss', 'val_loss'}
    if not required.issubset(rows[0]):
        raise ValueError(f'Missing history columns: {sorted(required - set(rows[0]))}')
    numeric = {}
    for field in rows[0]:
        try:
            values = [float(row[field]) for row in rows]
        except (TypeError, ValueError):
            continue
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f'Non-finite values in {field}')
        numeric[field] = values
    best_index = min(range(len(rows)), key=lambda index: numeric['val_loss'][index])
    summary = {
        'status': 'TRAINING_FINISHED_MODEL_QC_PENDING',
        'epochs_completed': len(rows),
        'best_epoch_zero_based': int(float(rows[best_index]['epoch'])),
        'best_val_loss': numeric['val_loss'][best_index],
        'final_val_loss': numeric['val_loss'][-1],
        'initial_val_loss': numeric['val_loss'][0],
        'best_training_loss': min(numeric['loss']),
        'history_columns': list(rows[0]),
        'scope': 'Finite-loss and convergence summary only; biological bias QC is pending.'
    }
    Path(output_path).write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    run(args.history, args.output)
