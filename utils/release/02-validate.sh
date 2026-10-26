#!/usr/bin/env bash

set -euo pipefail
# shellcheck source=utils/release/_common.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/_common.sh"

[[ $# -eq 0 ]] || die "usage: $0"
for command in cmake ctest git pre-commit; do
    require_command "$command"
done
require_clean_worktree

VERSION=$(current_version)
readonly VERSION
validate_version "$VERSION"
COMMIT=$(git rev-parse HEAD)
readonly COMMIT

pre-commit run --all-files

for preset in gcc-debug clang-debug gcc-release clang-release; do
    cmake --preset "$preset"
    cmake --build "build/$preset"
    ctest --test-dir "build/$preset" --output-on-failure -E '^cpptest$'
done

for preset in clang-asan clang-ubsan; do
    cmake --preset "$preset"
    cmake --build "build/$preset" --parallel 1
    ctest --test-dir "build/$preset" --output-on-failure -E '^cpptest$'
done

cmake --preset clang-tidy
cmake --build build/clang-tidy

ENGINE=$(container_engine)
readonly ENGINE
"$ENGINE" build --tag "chargefw:$VERSION" .

require_clean_worktree
[[ $(git rev-parse HEAD) == "$COMMIT" ]] || die "HEAD changed during validation"
printf 'Validated ChargeFW %s at %s.\n' "$VERSION" "$COMMIT"
