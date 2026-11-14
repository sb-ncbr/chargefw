#!/usr/bin/env bash

set -euo pipefail
# shellcheck source=utils/release/_common.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/_common.sh"

[[ $# -eq 0 ]] || die "usage: $0"
for command in cmake ctest git pre-commit systemd-run; do
    require_command "$command"
done

# Every build and test runs in a systemd scope limited to half of the physical memory, or to
# CHARGEFW_VALIDATE_MEMORY_MAX, without swap: a step that reaches the limit is stopped there instead
# of exhausting the machine. Tests are included because some of them build the project again.
systemd-run --user --scope --quiet true > /dev/null 2>&1 ||
    die "the validation runs its builds in memory-limited scopes and needs a systemd user session"
half_memory=$(awk '/^MemTotal:/ { print int($2 / 2048) "M" }' /proc/meminfo)
readonly MEMORY_MAX=${CHARGEFW_VALIDATE_MEMORY_MAX:-$half_memory}

# CHARGEFW_VALIDATE_JOBS sets the build parallelism, by default the build tool's own. It is passed
# as CMAKE_BUILD_PARALLEL_LEVEL, so the builds that tests start inherit it too.
if [[ -n ${CHARGEFW_VALIDATE_JOBS:-} ]]; then
    [[ $CHARGEFW_VALIDATE_JOBS =~ ^[1-9][0-9]*$ ]] ||
        die "CHARGEFW_VALIDATE_JOBS must be a positive integer, not '$CHARGEFW_VALIDATE_JOBS'"
    export CMAKE_BUILD_PARALLEL_LEVEL=$CHARGEFW_VALIDATE_JOBS
fi
printf 'Build jobs: %s, memory limit: %s\n' "${CMAKE_BUILD_PARALLEL_LEVEL:-default}" "$MEMORY_MAX"

limited() {
    systemd-run --user --scope --quiet -p MemoryMax="$MEMORY_MAX" -p MemorySwapMax=0 -- "$@"
}

require_clean_worktree

VERSION=$(current_version)
readonly VERSION
validate_version "$VERSION"
COMMIT=$(git rev-parse HEAD)
readonly COMMIT

pre-commit run --all-files

for preset in gcc-debug clang-debug gcc-release clang-release clang-asan clang-ubsan; do
    cmake --preset "$preset"
    limited cmake --build "build/$preset"
    limited ctest --test-dir "build/$preset" --output-on-failure -E '^cpptest$'
done

cmake --preset clang-tidy
limited cmake --build build/clang-tidy

ENGINE=$(container_engine)
readonly ENGINE
"$ENGINE" build --tag "chargefw:$VERSION" .

require_clean_worktree
[[ $(git rev-parse HEAD) == "$COMMIT" ]] || die "HEAD changed during validation"
printf 'Validated ChargeFW %s at %s.\n' "$VERSION" "$COMMIT"
