#!/usr/bin/env bash

set -euo pipefail
# shellcheck source=utils/release/_common.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/_common.sh"

[[ $# -eq 0 ]] || die "usage: $0"
for command in python3 sha256sum uv; do
    require_command "$command"
done
load_release_state
require_marker built
require_clean_worktree
[[ ! -f "$(marker_path local-verified)" ]] ||
    die "release artifacts are already verified; do not replace their checksum manifest"

shopt -s nullglob
artifacts=("$ARTIFACT_DIRECTORY"/*)
sdists=("$ARTIFACT_DIRECTORY"/chargefw-*.tar.gz)
wheels=("$ARTIFACT_DIRECTORY"/chargefw-*.whl)
[[ ${#artifacts[@]} -eq 5 ]] || die "expected exactly five release artifacts"
[[ ${#sdists[@]} -eq 1 ]] || die "expected exactly one source distribution"
[[ ${#wheels[@]} -eq 4 ]] || die "expected exactly four wheels"

for python_tag in 311 312 313 314; do
    matches=("$ARTIFACT_DIRECTORY"/chargefw-"$RELEASE_VERSION"-cp"$python_tag"-cp"$python_tag"-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl)
    [[ ${#matches[@]} -eq 1 ]] || die "expected one CPython $python_tag manylinux wheel"
done

for wheel in "${wheels[@]}"; do
    python3 -I -c 'import sys, zipfile; sys.exit(sys.argv[2] not in zipfile.ZipFile(sys.argv[1]).namelist())' \
        "$wheel" "chargefw-$RELEASE_VERSION.dist-info/licenses/licenses/README.md" ||
        die "$(basename -- "$wheel") does not contain the third-party license notices"
done

uvx twine check --strict "${artifacts[@]}"
(cd "$ARTIFACT_DIRECTORY" && LC_ALL=C sha256sum ./* | LC_ALL=C sort) > "$CHECKSUM_FILE"

readonly ENVIRONMENT="$RELEASE_ROOT/local-smoke"
uv venv --clear "$ENVIRONMENT" --python "${CHARGEFW_RELEASE_PYTHON:-3.14}"
uv pip install --python "$ENVIRONMENT/bin/python" 'numpy>=1.26'
uv pip install \
    --python "$ENVIRONMENT/bin/python" \
    --no-deps \
    --no-cache \
    --no-index \
    --find-links "$ARTIFACT_DIRECTORY" \
    --only-binary chargefw \
    "chargefw==$RELEASE_VERSION"
"$ENVIRONMENT/bin/python" - "$RELEASE_VERSION" <<'PY'
import sys
import chargefw

if chargefw.__version__ != sys.argv[1]:
    raise SystemExit(f"unexpected installed version: {chargefw.__version__}")
PY
"$ENVIRONMENT/bin/python" docs/recipes/calculate_file.py \
    tests/fixtures/synthetic/sdf/water.sdf \
    --format sdf \
    --method qeq \
    --parameter-set QEq_original

write_marker local-verified
printf 'Verified local artifacts and wrote %s.\n' "$CHECKSUM_FILE"
