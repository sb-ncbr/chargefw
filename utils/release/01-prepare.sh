#!/usr/bin/env bash

set -euo pipefail
# shellcheck source=utils/release/_common.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/_common.sh"

[[ $# -eq 1 ]] || die "usage: $0 VERSION"
readonly NEW_VERSION=$1
validate_version "$NEW_VERSION"
require_command git
require_command sed
require_clean_worktree

OLD_VERSION=$(current_version)
readonly OLD_VERSION
validate_version "$OLD_VERSION"
[[ "$NEW_VERSION" != "$OLD_VERSION" ]] || die "CMakeLists.txt already contains version $NEW_VERSION"
if git rev-parse --verify --quiet "refs/tags/v$NEW_VERSION" >/dev/null; then
    die "tag v$NEW_VERSION already exists"
fi

[[ "$NEW_VERSION" =~ ^([0-9]+\.[0-9]+\.[0-9]+)(.*)$ ]]
readonly NEW_NUMERIC_VERSION=${BASH_REMATCH[1]}
readonly NEW_VERSION_SUFFIX=${BASH_REMATCH[2]}
sed -i -E \
    -e "s/^([[:space:]]*VERSION[[:space:]]+)[0-9]+\.[0-9]+\.[0-9]+([[:space:]]*)$/\1${NEW_NUMERIC_VERSION}\2/" \
    -e "s/^set\(CHARGEFW_VERSION_SUFFIX \"[^\"]*\"\)$/set(CHARGEFW_VERSION_SUFFIX \"${NEW_VERSION_SUFFIX}\")/" \
    CMakeLists.txt

[[ $(current_version) == "$NEW_VERSION" ]] || die "version update failed"
printf 'Updated ChargeFW %s -> %s. Review and commit CMakeLists.txt, then run 02-validate.sh.\n' \
    "$OLD_VERSION" "$NEW_VERSION"
