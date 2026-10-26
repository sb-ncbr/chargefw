#!/usr/bin/env bash

set -euo pipefail

REPOSITORY_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
readonly REPOSITORY_ROOT

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

container_engine() {
    if [[ -n ${CHARGEFW_CONTAINER_ENGINE:-} ]]; then
        printf '%s\n' "$CHARGEFW_CONTAINER_ENGINE"
    elif command -v podman >/dev/null 2>&1; then
        printf 'podman\n'
    elif command -v docker >/dev/null 2>&1; then
        printf 'docker\n'
    else
        die "Podman or Docker is required"
    fi
}
