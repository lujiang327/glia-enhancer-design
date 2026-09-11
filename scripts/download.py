#!/usr/bin/env python3
"""Fetch manifest-selected public assets, checksum, then atomically rename.

Usage: python3 scripts/download.py config/assets.json
Existing valid files are skipped. Partial files are never accepted as complete.
"""
import hashlib
import json
from pathlib import Path
import sys
import urllib.request


def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(2**20),b''):h.update(chunk)
    return h.hexdigest()


def fetch(asset):
    path=Path(asset['path']);path.parent.mkdir(parents=True,exist_ok=True)
    if path.is_file() and path.stat().st_size==asset['bytes'] and sha256(path)==asset['sha256']:
        print('Verified existing:',path);return
    partial=path.with_suffix(path.suffix+'.part')
    with urllib.request.urlopen(asset['url'],timeout=60) as source, partial.open('wb') as target:
        while True:
            block=source.read(2**20)
            if not block:break
            target.write(block)
    if partial.stat().st_size!=asset['bytes'] or sha256(partial)!=asset['sha256']:
        raise ValueError(f'Checksum or size failure: {partial}')
    partial.replace(path);print('Downloaded and verified:',path)


if __name__=='__main__':
    for asset in json.loads(Path(sys.argv[1]).read_text()):fetch(asset)
