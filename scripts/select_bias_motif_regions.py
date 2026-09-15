#!/usr/bin/env python3
"""Select a deterministic, chromosome-held-out region subset for bias motifs."""
import argparse
import hashlib
import json
from pathlib import Path
import random


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def run(regions_path, fold_path, split_name, number, seed, output_path, manifest_path):
    fold = json.loads(Path(fold_path).read_text())
    chromosomes = set(fold[split_name])
    with open(regions_path) as handle:
        eligible = [line for line in handle if line.split('\t', 1)[0] in chromosomes]
    if len(eligible) < number:
        raise ValueError(f'Requested {number} regions from only {len(eligible)} eligible regions')

    chosen = sorted(random.Random(seed).sample(range(len(eligible)), number))
    selected = [eligible[index] for index in chosen]
    Path(output_path).write_text(''.join(selected))
    selected_chromosomes = sorted({line.split('\t', 1)[0] for line in selected})
    manifest = {
        'status': 'PASS',
        'source': str(regions_path),
        'fold': str(fold_path),
        'split': split_name,
        'eligible_regions': len(eligible),
        'selected_regions': len(selected),
        'seed': seed,
        'chromosomes': selected_chromosomes,
        'output': str(output_path),
        'output_sha256': sha256(output_path),
    }
    Path(manifest_path).write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--regions', required=True)
    parser.add_argument('--fold', required=True)
    parser.add_argument('--split', default='test')
    parser.add_argument('--number', type=int, default=30000)
    parser.add_argument('--seed', type=int, default=1234)
    parser.add_argument('--output', required=True)
    parser.add_argument('--manifest', required=True)
    args = parser.parse_args()
    run(args.regions, args.fold, args.split, args.number, args.seed, args.output, args.manifest)
