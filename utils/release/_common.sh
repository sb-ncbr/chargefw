#!/usr/bin/env bash

set -euo pipefail

RELEASE_SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
readonly RELEASE_SCRIPT_DIR
REPOSITORY_ROOT=$(cd -- "$RELEASE_SCRIPT_DIR/../.." && pwd)
readonly REPOSITORY_ROOT
readonly RELEASE_ROOT="$REPOSITORY_ROOT/build/release"
readonly RELEASE_STATE_FILE="$RELEASE_ROOT/state.env"

cd "$REPOSITORY_ROOT"

die() {
    printf 'error: %s\n' "$*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

current_version() {
    local numeric suffix
    numeric=$(sed -nE 's/^[[:space:]]*VERSION[[:space:]]+([0-9]+\.[0-9]+\.[0-9]+)[[:space:]]*$/\1/p' \
        CMakeLists.txt)
    suffix=$(sed -nE 's/^set\(CHARGEFW_VERSION_SUFFIX "([^"]*)"\)$/\1/p' CMakeLists.txt)
    printf '%s%s\n' "$numeric" "$suffix"
}

validate_version() {
    # Canonical PEP 440 spelling keeps artifact filenames identical to the recorded version.
    local number='(0|[1-9][0-9]*)'
    [[ "$1" =~ ^$number\.$number\.$number((a|b|rc)$number)?$ ]] ||
        die "version must have the form MAJOR.MINOR.PATCH with an optional aN, bN, or rcN suffix" \
            "and no leading zeros"
}

require_clean_worktree() {
    [[ -z $(git status --porcelain --untracked-files=all) ]] ||
        die "the Git worktree must be clean"
}

initialize_release_state() {
    local version=$1
    local commit=$2

    rm -rf -- "$RELEASE_ROOT"
    mkdir -p -- "$RELEASE_ROOT"
    {
        printf 'RELEASE_VERSION=%q\n' "$version"
        printf 'RELEASE_COMMIT=%q\n' "$commit"
        printf 'ARTIFACT_DIRECTORY=%q\n' "build/release/artifacts"
        printf 'CHECKSUM_FILE=%q\n' "build/release/sha256sums.txt"
    } > "$RELEASE_STATE_FILE"
}

load_release_state() {
    [[ -f "$RELEASE_STATE_FILE" ]] ||
        die "release state is missing; run 02-validate.sh first"

    # This file is generated locally by initialize_release_state.
    # shellcheck disable=SC1090
    source "$RELEASE_STATE_FILE"
    readonly RELEASE_VERSION RELEASE_COMMIT ARTIFACT_DIRECTORY CHECKSUM_FILE

    validate_version "$RELEASE_VERSION"
    [[ $(current_version) == "$RELEASE_VERSION" ]] ||
        die "CMakeLists.txt no longer contains release version $RELEASE_VERSION"
    [[ $(git rev-parse HEAD) == "$RELEASE_COMMIT" ]] ||
        die "HEAD no longer points to release commit $RELEASE_COMMIT"
}

marker_path() {
    printf '%s/%s.ok\n' "$RELEASE_ROOT" "$1"
}

write_marker() {
    local checksum_digest=
    if [[ -f ${CHECKSUM_FILE:-} ]]; then
        checksum_digest=$(sha256sum "$CHECKSUM_FILE")
        checksum_digest=${checksum_digest%% *}
    fi
    printf '%s %s %s\n' "$RELEASE_VERSION" "$RELEASE_COMMIT" "$checksum_digest" \
        > "$(marker_path "$1")"
}

require_marker() {
    local marker
    local marker_version marker_commit marker_checksum
    marker=$(marker_path "$1")
    [[ -f "$marker" ]] || die "release gate '$1' has not completed"
    read -r marker_version marker_commit marker_checksum < "$marker"
    [[ "$marker_version" == "$RELEASE_VERSION" && "$marker_commit" == "$RELEASE_COMMIT" ]] ||
        die "release gate '$1' belongs to another release"
    if [[ -n "$marker_checksum" ]]; then
        [[ -f "$CHECKSUM_FILE" ]] || die "release checksum manifest is missing"
        local checksum_digest
        checksum_digest=$(sha256sum "$CHECKSUM_FILE")
        checksum_digest=${checksum_digest%% *}
        [[ "$checksum_digest" == "$marker_checksum" ]] ||
            die "release checksum manifest changed after gate '$1'"
    fi
}

confirm_release_action() {
    local action=$1
    local confirmation=${CHARGEFW_RELEASE_CONFIRM:-}

    if [[ "$confirmation" != "$RELEASE_VERSION" ]]; then
        [[ -t 0 ]] || die "set CHARGEFW_RELEASE_CONFIRM=$RELEASE_VERSION to confirm $action"
        printf '%s for ChargeFW %s. Type the version to continue: ' "$action" "$RELEASE_VERSION" >&2
        read -r confirmation
    fi
    [[ "$confirmation" == "$RELEASE_VERSION" ]] || die "$action was not confirmed"
}

require_artifacts() {
    require_command sha256sum
    [[ -d "$ARTIFACT_DIRECTORY" ]] || die "artifact directory is missing"
    [[ -f "$CHECKSUM_FILE" ]] || die "artifact checksum manifest is missing"
    (cd "$ARTIFACT_DIRECTORY" && sha256sum --check "$REPOSITORY_ROOT/$CHECKSUM_FILE")
}

container_engine() {
    if [[ -n ${CIBW_CONTAINER_ENGINE:-} ]]; then
        printf '%s\n' "$CIBW_CONTAINER_ENGINE"
    elif command -v podman >/dev/null 2>&1; then
        printf 'podman\n'
    elif command -v docker >/dev/null 2>&1; then
        printf 'docker\n'
    else
        die "Podman or Docker is required"
    fi
}

remote_artifact_statuses() {
    local json_url=$1
    local output_file=$2

    require_command python3
    python3 - "$json_url" "$CHECKSUM_FILE" "$ARTIFACT_DIRECTORY" "$output_file" <<'PY'
import json
from pathlib import Path
import sys
from urllib.error import HTTPError
from urllib.request import Request, urlopen

json_url, checksum_path, artifact_directory, output_path = sys.argv[1:]
request = Request(json_url, headers={"User-Agent": "chargefw-release-script"})
try:
    with urlopen(request, timeout=30) as response:
        payload = json.load(response)
except HTTPError as error:
    if error.code != 404:
        raise
    payload = {"urls": []}

remote = {entry["filename"]: entry["digests"]["sha256"] for entry in payload["urls"]}
statuses = []
failed = False
local_filenames = set()
for line in Path(checksum_path).read_text(encoding="utf-8").splitlines():
    digest, artifact_path = line.split(maxsplit=1)
    filename = Path(artifact_path.lstrip("*")).name
    artifact_path = str(Path(artifact_directory) / filename)
    local_filenames.add(filename)
    remote_digest = remote.get(filename)
    if remote_digest is None:
        status = "missing"
    elif remote_digest == digest:
        status = "match"
    else:
        status = "mismatch"
        failed = True
    statuses.append(f"{status}\t{artifact_path}")

unexpected = sorted(remote.keys() - local_filenames)
if unexpected:
    print("remote index contains unexpected artifacts: " + ", ".join(unexpected), file=sys.stderr)
    failed = True

Path(output_path).write_text("\n".join(statuses) + "\n", encoding="utf-8")
if failed:
    raise SystemExit("a remote artifact has the expected filename but a different SHA-256")
PY
}

upload_missing_artifacts() {
    local repository=$1
    local repository_url=$2
    local json_url=$3
    local status_file="$RELEASE_ROOT/upload-status.txt"
    local -a missing=()
    local status artifact

    remote_artifact_statuses "$json_url" "$status_file"
    while IFS=$'\t' read -r status artifact; do
        case "$status" in
            missing) missing+=("$artifact") ;;
            match) printf 'Already uploaded with matching hash: %s\n' "$(basename -- "$artifact")" ;;
            *) die "unexpected remote artifact status: $status" ;;
        esac
    done < "$status_file"

    if ((${#missing[@]} == 0)); then
        printf 'All artifacts are already present with matching hashes.\n'
        return
    fi

    local -a command=(env -u TWINE_REPOSITORY -u TWINE_REPOSITORY_URL
                      TWINE_USERNAME=__token__ uvx twine upload
                      --repository "$repository" --repository-url "$repository_url")
    command+=("${missing[@]}")
    "${command[@]}"
}

preflight_remote_artifacts() {
    local json_url=$1
    local status_file="$RELEASE_ROOT/preflight-status.txt"

    remote_artifact_statuses "$json_url" "$status_file"
    printf 'Remote artifact namespace has no conflicting files.\n'
}

verify_remote_artifacts() {
    local json_url=$1
    local status_file="$RELEASE_ROOT/verify-status.txt"
    local status artifact

    remote_artifact_statuses "$json_url" "$status_file"
    while IFS=$'\t' read -r status artifact; do
        [[ "$status" == match ]] ||
            die "remote artifact is not available with the expected hash: $(basename -- "$artifact")"
    done < "$status_file"
}

smoke_test_index() {
    local json_url=$1
    local environment=$2
    local python=${CHARGEFW_RELEASE_PYTHON:-3.14}
    local downloaded_wheel

    uv venv --clear "$environment" --python "$python"
    uv pip install --python "$environment/bin/python" 'numpy>=1.26'
    downloaded_wheel=$("$environment/bin/python" - "$json_url" "$CHECKSUM_FILE" \
        "$environment" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import Request, urlopen

json_url, checksum_path, environment = sys.argv[1:]
python_tag = f"cp{sys.version_info.major}{sys.version_info.minor}"
expected = {}
for line in Path(checksum_path).read_text(encoding="utf-8").splitlines():
    digest, artifact_path = line.split(maxsplit=1)
    expected[Path(artifact_path.lstrip("*")).name] = digest

request = Request(json_url, headers={"User-Agent": "chargefw-release-script"})
with urlopen(request, timeout=30) as response:
    payload = json.load(response)
matches = [
    entry
    for entry in payload["urls"]
    if entry["filename"].endswith(".whl") and f"-{python_tag}-{python_tag}-" in entry["filename"]
]
if len(matches) != 1:
    raise SystemExit(f"expected one remote {python_tag} wheel, found {len(matches)}")
entry = matches[0]
if expected.get(entry["filename"]) != entry["digests"]["sha256"]:
    raise SystemExit("remote wheel is not present in the approved checksum manifest")

wheel_path = Path(environment) / entry["filename"]
with urlopen(Request(entry["url"], headers={"User-Agent": "chargefw-release-script"}), timeout=60) as response:
    wheel = response.read()
if hashlib.sha256(wheel).hexdigest() != entry["digests"]["sha256"]:
    raise SystemExit("downloaded wheel SHA-256 does not match the package index")
wheel_path.write_bytes(wheel)
print(wheel_path)
PY
    )
    uv pip install \
        --python "$environment/bin/python" \
        --no-deps \
        --no-cache \
        --no-index \
        "$downloaded_wheel"

    "$environment/bin/python" - "$RELEASE_VERSION" <<'PY'
import sys
import chargefw

expected = sys.argv[1]
if chargefw.__version__ != expected:
    raise SystemExit(f"expected chargefw {expected}, imported {chargefw.__version__}")
PY
    "$environment/bin/python" docs/recipes/calculate_file.py \
        tests/fixtures/synthetic/sdf/water.sdf \
        --format sdf \
        --method qeq \
        --parameter-set QEq_original
}
