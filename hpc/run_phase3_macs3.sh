#!/usr/bin/env bash
set -euo pipefail

# Run the scratch-installed, pinned MACS3 environment inside the checksum-pinned
# ChromBPNet container. PHASE3_SCRATCH_ROOT must be supplied by the caller.

if [[ "$#" -eq 0 ]]; then
  echo "Usage: PHASE3_SCRATCH_ROOT=... $0 COMMAND [ARG ...]" >&2
  exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$script_dir/.." && pwd)
scratch_root="${PHASE3_SCRATCH_ROOT:?Set PHASE3_SCRATCH_ROOT}"
image="$project_root/containers/chrombpnet-upstream-20260914.sif"
environment="$scratch_root/software/macs3-3.0.4"
expected_image_sha256="3eea58b1606fe1245ad70abe70ea6268cf10587ea88de52263c268d1952249d5"

[[ -f "$image" ]] || { echo "Missing container: $image" >&2; exit 1; }
[[ -x "$environment/bin/python" ]] || { echo "Missing MACS3 environment: $environment" >&2; exit 1; }
command -v singularity >/dev/null || { echo "Load the Singularity module first" >&2; exit 1; }

observed_image_sha256=$(sha256sum "$image" | awk '{print $1}')
[[ "$observed_image_sha256" == "$expected_image_sha256" ]] || {
  echo "Container SHA256 mismatch: $observed_image_sha256" >&2
  exit 1
}

exec singularity exec \
  --cleanenv \
  --no-mount hostfs \
  --bind "$project_root:/workspace" \
  --bind "$scratch_root:/phase3" \
  "$image" \
  bash -c '
    set -euo pipefail
    export PATH=/phase3/software/macs3-3.0.4/bin:/opt/conda/bin:/usr/local/bin:/usr/bin:/bin
    cd /workspace
    exec "$@"
  ' bash "$@"
