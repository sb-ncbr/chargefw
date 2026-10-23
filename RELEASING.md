# ChargeFW Maintainer Release Procedure

This procedure validates the repository, builds immutable Python artifacts, promotes them through
TestPyPI, tags their source commit, and publishes the same files to PyPI. The numbered scripts under
`utils/release/` enforce that order and retain release state under the ignored `build/release/` directory.

## Current wheel matrix

The current release target is CPython 3.11 through 3.14 on Linux x86-64. `cibuildwheel` builds each wheel
in a manylinux container, repairs its native libraries, installs it, and runs the Python test suite. The
resulting wheels carry `manylinux_2_27_x86_64` and `manylinux_2_28_x86_64` tags.

CPython 3.15 is not included because the pinned optional Gemmi 0.7.4 package does not currently provide a
usable CPython 3.15 wheel. Other operating systems and architectures are not yet qualified.

## Requirements and safety

Install the development dependencies from [DEVELOPMENT.md](DEVELOPMENT.md), `sed`, Python 3, `uv`,
`pre-commit`, and either Podman or Docker. Set `CIBW_CONTAINER_ENGINE=docker` to choose Docker when both
engines are installed; otherwise the scripts prefer Podman.
The optional `keyring` command can retrieve stored upload tokens without typing them at each release.

The release version is `MAJOR.MINOR.PATCH` with an optional PEP 440 pre-release suffix `aN`, `bN`, or
`rcN`, for example `3.0.0a1`. In `CMakeLists.txt`, the numeric part is the `project(... VERSION ...)`
declaration and the suffix is `CHARGEFW_VERSION_SUFFIX`, empty for final releases. The full version is
the Python package version, `chargefw.__version__`, the generator version in output files, and the tag
name; the CMake package version is the numeric part. `01-prepare.sh` updates both parts.

Each script verifies the recorded version and Git commit. Publication scripts additionally verify the
exact artifact checksums. A repeated upload skips an existing remote file only when its SHA-256 matches;
the release stops if a filename exists with different content.

Published files cannot be replaced. If any artifact must be rebuilt after it has been uploaded to either
package index, increment the version and restart the release. Do not delete or reuse a published version
or move a published release tag.

## 1. Prepare the release commit

Start from an up-to-date branch and a clean worktree. Select the version explicitly:

```bash
git fetch origin
git switch -c release/0.1.3 origin/main
./utils/release/01-prepare.sh 0.1.3
git diff -- CMakeLists.txt
git add CMakeLists.txt
git commit -m "Prepare release 0.1.3"
```

Include any release notes or other intentional release changes in that commit. The remaining phases
require a clean worktree.

## 2. Validate and build

Run the complete GCC, Clang, release, sanitizer, and clang-tidy matrix and build the distributable CLI
container image:

```bash
./utils/release/02-validate.sh
```

The sanitizer builds run sequentially with one build job. The full preset tests include native, CLI,
Python, installed-package, relocation, and downstream CMake consumer checks.

Build one source distribution and then build all wheels from that source distribution:

```bash
./utils/release/03-build.sh
./utils/release/04-verify-local.sh
```

The local verification requires exactly one sdist and the four configured wheels, runs strict package
metadata checks, writes `build/release/sha256sums.txt`, installs a wheel in a clean environment, verifies
its runtime version, and runs a smoke calculation. Do not rebuild after this point.

## 3. Publish the source commit

Fast-forward `main` to the validated release commit, then let the gated script push it:

```bash
git switch main
git merge --ff-only release/0.1.3
./utils/release/05-push-main.sh
```

The script requires the current branch to be `main`, asks for the release version as confirmation, pushes
only `main`, and verifies that `origin/main` points to the recorded release commit.

## 4. TestPyPI

Use a TestPyPI project-scoped token. To store it in the system keyring, run `keyring set` once; it prompts
for the token without putting it in shell history:

```bash
keyring set 'https://test.pypi.org/legacy/' __token__
```

`uvx twine` runs in an isolated environment that may not see the system's keyring backend (for example,
KWallet). Read the token through the system `keyring` command and pass it to the upload script without
displaying it or placing it in shell history:

```bash
TWINE_PASSWORD="$(keyring get 'https://test.pypi.org/legacy/' __token__)" \
  ./utils/release/06-publish-testpypi.sh
```

Do not run the command with shell tracing (`set -x`) enabled. If no keyring is available, use a masked
prompt instead:

```bash
read -rsp 'TestPyPI token: ' TWINE_PASSWORD && export TWINE_PASSWORD
printf '\n'
./utils/release/06-publish-testpypi.sh
unset TWINE_PASSWORD
```

Verify every remote hash, install the exact wheel from TestPyPI in a clean environment, assert the
runtime version, and run the smoke calculation:

```bash
./utils/release/07-verify-testpypi.sh
```

Package-index propagation can take time. If verification reports a missing artifact, wait and rerun only
`07-verify-testpypi.sh`; do not rebuild or upload a second set of files.

## 5. Tag the approved source

Create and push the annotated tag only after TestPyPI verification:

```bash
./utils/release/08-tag.sh
```

The script creates `vVERSION` at the recorded release commit, pushes it, and verifies that the remote tag
resolves to that commit. Before creating the tag it also checks production PyPI for conflicting artifact
filenames. If a local tag already exists, it must be annotated and point to the same commit.

## 6. Production PyPI

Use a separate production PyPI project token. Store it under the production URL if needed:

```bash
keyring set 'https://upload.pypi.org/legacy/' __token__
```

Retrieve it from the system keyring for the production upload. The production step consumes the unchanged
local checksum manifest approved by TestPyPI:

```bash
TWINE_PASSWORD="$(keyring get 'https://upload.pypi.org/legacy/' __token__)" \
  ./utils/release/09-publish-pypi.sh
```

The same masked-prompt fallback applies when a system keyring is unavailable:

```bash
read -rsp 'PyPI token: ' TWINE_PASSWORD && export TWINE_PASSWORD
printf '\n'
./utils/release/09-publish-pypi.sh
unset TWINE_PASSWORD
```

Finally, compare every PyPI hash and repeat the clean installation and smoke calculation:

```bash
./utils/release/10-verify-pypi.sh
```

If PyPI propagation is incomplete, rerun only `10-verify-pypi.sh`. A successful run writes the
`build/release/release-complete.ok` gate.

## Recovery and retention

The scripts are resumable. Rerun a failed phase after correcting an environmental or network problem.
Completed build and local-verification phases deliberately refuse to run again so they cannot replace an
approved artifact set. Upload phases query the target index before uploading, which permits safe recovery
from a partial upload when every file already present has the expected hash.

If source code, the version, the release commit, or an artifact changes, restart from the appropriate
earlier phase. Once an upload has occurred, a required rebuild means creating a new version.

Keep release artifacts and records in `build/release/`. Do not create optional backup archives or copy
release files into the repository's parent directory. The next release validation replaces state for a
different version only when the previous release has a `release-complete.ok` gate; no backup step is
required before replacement.

This output-location rule applies to project build outputs and release records, not tool-managed caches,
container storage, or runtime files. Let Podman and other build tools use their normal configured
locations, including the system-managed `XDG_RUNTIME_DIR`; do not relocate those directories into the
repository to satisfy this rule.

Create the GitHub release from the pushed tag and attach the release notes. Release automation, trusted
publishing, durable CI artifact retention, and additional wheel platforms remain tracked in
[TODO.md](TODO.md#distribution).
