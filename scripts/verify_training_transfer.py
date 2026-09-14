#!/usr/bin/env python3
"""Verify transferred input sizes and SHA256 values using Python stdlib only."""
import hashlib
import json
from pathlib import Path


def run():
    records = json.loads(Path('reports/preparation/training_input_manifest.json').read_text())
    for record in records:
        path = Path(record['path'])
        if path.stat().st_size != record['bytes']:
            raise ValueError(f'Size mismatch: {path}')
        checksum = hashlib.sha256()
        with path.open('rb') as handle:
            for block in iter(lambda: handle.read(2**20), b''):
                checksum.update(block)
        if checksum.hexdigest() != record['sha256']:
            raise ValueError(f'SHA256 mismatch: {path}')
        print(f'PASS {path}', flush=True)
    print(f'Verified {len(records)} training inputs. Model validation remains separate.')


if __name__ == '__main__':
    run()
