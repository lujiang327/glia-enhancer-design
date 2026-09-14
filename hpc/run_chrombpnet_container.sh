#!/usr/bin/env bash
set -euo pipefail

# Run a command in the checksum-pinned container with the checksum-pinned
# ChromBPNet source taking precedence over the package bundled in the image.
# Usage: run_chrombpnet_container.sh [--gpu] COMMAND [ARG ...]

use_gpu=0
if [[ "${1:-}" == "--gpu" ]]; then
  use_gpu=1
  shift
fi
if [[ "$#" -eq 0 ]]; then
  echo "Usage: $0 [--gpu] COMMAND [ARG ...]" >&2
  exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$script_dir/.." && pwd)
image="$project_root/containers/chrombpnet-upstream-20260914.sif"
source_dir="$project_root/data/intermediate/vendor/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c"
expected_image_sha256="3eea58b1606fe1245ad70abe70ea6268cf10587ea88de52263c268d1952249d5"

[[ -f "$image" ]] || { echo "Missing container: $image" >&2; exit 1; }
[[ -f "$source_dir/setup.py" ]] || { echo "Missing pinned source: $source_dir" >&2; exit 1; }
command -v singularity >/dev/null || { echo "Load the Singularity module first" >&2; exit 1; }

observed_image_sha256=$(sha256sum "$image" | awk '{print $1}')
if [[ "$observed_image_sha256" != "$expected_image_sha256" ]]; then
  echo "Container SHA256 mismatch: $observed_image_sha256" >&2
  exit 1
fi

singularity_args=(
  exec
  --cleanenv
  --no-mount hostfs
  --bind "$project_root:/workspace"
)
if [[ "$use_gpu" -eq 1 ]]; then
  singularity_args+=(--nv)
fi

cuda_visible_devices="${CUDA_VISIBLE_DEVICES:-}"
exec singularity "${singularity_args[@]}" "$image" \
  bash -c '
    set -euo pipefail
    export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin
    export PYTHONPATH=/workspace/data/intermediate/vendor/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c
    export LD_LIBRARY_PATH=/.singularity.d/libs:/usr/local/cuda/lib64:/usr/local/cuda/extras/CUPTI/lib64
    export CUDA_VISIBLE_DEVICES="$1"
    shift
    cd /workspace
    exec "$@"
  ' bash "$cuda_visible_devices" "$@"
