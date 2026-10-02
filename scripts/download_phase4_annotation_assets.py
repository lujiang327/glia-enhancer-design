#!/usr/bin/env python3
"""Snapshot UCSC annotation assets and refuse changes to a saved snapshot."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request


def digest(path):
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    root = Path(config["assets_directory"])
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "asset_snapshot.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {"assets": {}}
    for asset in config["assets"]:
        name, url = asset["name"], asset["url"]
        path = root / name
        if name in manifest["assets"]:
            saved = manifest["assets"][name]
            if saved["url"] != url or not path.exists() or digest(path) != saved["sha256"]:
                raise ValueError("Saved asset snapshot mismatch: {}".format(name))
            print("Verified {}".format(name), flush=True)
            continue
        if path.exists():
            raise ValueError("Unmanifested asset exists; review before reuse: {}".format(path))
        partial = path.with_name(path.name + ".part")
        print("Downloading {}".format(url), flush=True)
        hasher = hashlib.sha256()
        total = 0
        with urllib.request.urlopen(url, timeout=120) as response, partial.open("wb") as target:
            last_modified = response.headers.get("Last-Modified")
            etag = response.headers.get("ETag")
            expected_length = response.headers.get("Content-Length")
            resolved_url = response.geturl()
            for chunk in iter(lambda: response.read(1024 * 1024), b""):
                target.write(chunk)
                hasher.update(chunk)
                total += len(chunk)
        if total == 0 or (expected_length is not None and total != int(expected_length)):
            raise ValueError("Incomplete annotation asset: {}".format(name))
        partial.replace(path)
        manifest["assets"][name] = {
            "url": url, "resolved_url": resolved_url, "sha256": hasher.hexdigest(),
            "bytes": total, "http_last_modified": last_modified, "http_etag": etag,
            "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
            "verification_scope": "observed snapshot checksum, not an upstream published checksum",
        }
        manifest_partial = manifest_path.with_suffix(".json.part")
        manifest_partial.write_text(json.dumps(manifest, indent=2) + "\n")
        manifest_partial.replace(manifest_path)
        print("Saved {} bytes SHA256 {}".format(total, hasher.hexdigest()), flush=True)


if __name__ == "__main__":
    main()
