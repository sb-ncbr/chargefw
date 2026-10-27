#!/usr/bin/env bash
# Installs one release wheel in a clean environment, checks its version, and runs a calculation.

set -euo pipefail

[[ $# -eq 2 ]] || { printf 'usage: %s VERSION WHEEL\n' "$0" >&2; exit 2; }
version=$1
wheel=$(realpath -- "$2")
repository=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
environment=$(mktemp -d)
trap 'rm -rf -- "$environment"' EXIT

python3 -m venv "$environment"
"$environment/bin/python" -m pip install --quiet 'numpy>=2.0'
"$environment/bin/python" -m pip install --quiet --no-deps --no-index "$wheel"
# Run outside the checkout so that only the installed package can be imported.
cd "$environment"
"$environment/bin/python" -c \
    'import sys, chargefw; sys.exit(chargefw.__version__ != sys.argv[1])' "$version"
"$environment/bin/python" "$repository/docs/examples/calculate_file.py" \
    "$repository/tests/fixtures/synthetic/sdf/water.sdf" \
    --format sdf --method qeq --parameter-set QEq_original
