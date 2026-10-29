# ChargeFW Maintainer Release Procedure

A release is a reviewed version commit on `main`, validated locally, and published by the
[release workflow](.github/workflows/release.yml). The workflow builds the artifacts once, publishes
them to TestPyPI, tags the commit, and publishes the same files to PyPI after maintainer approval. It
then creates the GitHub Release and pushes the CLI image to Docker Hub.

## Current wheel matrix

The release publishes the sdist and one wheel per platform. Each wheel is built against CPython's
stable ABI (`cp312-abi3`), so it serves CPython 3.12 and later.

| Platform | Wheel platform tag |
| --- | --- |
| Linux x86-64 | `manylinux_2_27_x86_64.manylinux_2_28_x86_64` |
| Linux aarch64 | `manylinux_2_26_aarch64.manylinux_2_28_aarch64` |
| macOS 14 or newer, Apple silicon | `macosx_14_0_arm64` |
| Windows x86-64 | `win_amd64` |

`.github/workflows/wheels.yml`, which CI also runs, builds the sdist once and each wheel from it on a
runner of that platform with `cibuildwheel`. It repairs the native libraries, and installs the wheel and
runs the Python test suite with CPython 3.12, 3.13, and 3.14. The build settings are in `pyproject.toml`.
auditwheel derives the manylinux versions from the glibc functions the binaries use, so they can change
between releases; the release workflow checks for one wheel per platform, whatever its manylinux tag.

The wheels also install on CPython 3.15, but they are not tested there, because the pinned optional Gemmi
0.7.4 package does not currently provide a usable CPython 3.15 wheel.

## Versions

The release version is `MAJOR.MINOR.PATCH` with an optional PEP 440 pre-release suffix `aN`, `bN`, or
`rcN`, for example `3.0.0a2`. In `CMakeLists.txt`, the numeric part is the `project(... VERSION ...)`
declaration and the suffix is `CHARGEFW_VERSION_SUFFIX`, empty for final releases. The full version is
the Python package version, `chargefw.__version__`, the generator version in output files, and the tag
name `vVERSION`; the CMake package version is the numeric part.

Published files cannot be replaced. If an artifact must change after it has been uploaded to either
package index, release a new version. Do not delete or reuse a published version or move a release tag.

## 1. Prepare the version commit

Start from an up-to-date `main` and a clean worktree:

```bash
./utils/release/01-prepare.sh 3.0.0a3
```

The script updates both version parts in `CMakeLists.txt` and the Docker image tag in `README.md` and
`docs/getting-started/installation.md`. Rename the `## Unreleased` heading in `CHANGELOG.md` to
`## 3.0.0a3`; the workflow requires that section and uses it as the GitHub Release notes, so its links
must be absolute. Commit these files, for example as "Prepare release 3.0.0a3".

## 2. Validate locally

CI tests `gcc-debug` and the wheels. Run the complete GCC, Clang, release, sanitizer, and clang-tidy
matrix and build the distributable CLI image on the version commit:

```bash
./utils/release/02-validate.sh
```

It requires `cmake`, `ctest`, `pre-commit`, the development dependencies from
[DEVELOPMENT.md](DEVELOPMENT.md), and Podman or Docker; set `CHARGEFW_CONTAINER_ENGINE=docker` to choose
Docker when both are installed. The sanitizer builds run sequentially with one build job.

The tests check generated charge categories against the charges dictionary. Also check generated mmCIF
against the current PDBx dictionary, with an installed `chargefw` and the Gemmi command-line program:

```bash
dir=$(mktemp -d)
curl -sSfL -o "$dir/mmcif_pdbx_v50.dic" https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic
chargefw calculate --method formal --output-mmcif tests/fixtures/corpus/sdf/audit/representative.sdf "$dir/sdf"
chargefw calculate --method formal --output-mmcif tests/fixtures/corpus/cif/1ek9.cif "$dir/cif"
gemmi validate -d "$dir/mmcif_pdbx_v50.dic" -d tests/fixtures/dictionaries/mmcif_charges_v11.dic \
  "$dir"/sdf/*.cif "$dir"/cif/*.cif
```

It should exit with status 0. Notes that `mmcif_charges_v11.dic` does not match `mmcif_charges.dic` are
expected: the charges dictionary is declared under its file name, as the charges-schema project
prescribes. Parent-category checks (`-p`) are not expected to pass, because the output omits categories
such as `_entity` and `_atom_type`.

## 3. Push and wait for CI

Push the version commit to `main` and wait until `ci.yml` succeeds on it. The release workflow refuses
to run on a commit without a successful CI run.

## 4. Run the release workflow

```bash
gh workflow run release.yml -f version=3.0.0a3
```

or use **Actions → Release → Run workflow** on `main`. The jobs run in this order:

| Job | What it does |
| --- | --- |
| `qualify` | Checks that the run is on `main`, the version matches `CMakeLists.txt` and has a changelog section, the tag and both index versions are unused, and CI passed |
| `wheels` | Builds the sdist and one wheel per platform from it, and tests each wheel, as CI does |
| `build` | Checks the artifact set, notices, and metadata, runs a smoke calculation, and stores the artifacts with a `SHA256SUMS` manifest for 90 days |
| `testpypi` | Uploads to TestPyPI, compares the remote hashes, and runs the smoke calculation with the downloaded wheel |
| `tag` | Checks PyPI for conflicting files, then pushes the annotated tag `vVERSION` |
| `pypi` | Waits for approval, then uploads the same files to PyPI and verifies them as for TestPyPI |
| `github-release` | Creates the GitHub Release with the artifacts, manifest, and changelog section |
| `docker-image` | Pushes `sbncbr/chargefw:VERSION`; final releases also move `MAJOR.MINOR` and `latest` |

When the run stops at `pypi`, check the TestPyPI project page, then approve the deployment in the run
(**Review deployments**).

## 5. After the release

- Check that the release appears as a new version of the Figshare item and that Read the Docs built the
  `vVERSION` documentation.
- Add a new `## Unreleased` section to `CHANGELOG.md` with the next user-visible change.

## Recovery

Each job can be resumed. After a network, index, or registry failure, use **Re-run failed jobs**; it
reuses the stored artifacts. Upload jobs skip a remote file only when its SHA-256 matches the manifest
and stop if a filename exists with different content. Package indexes can return an outdated file list
for a short time after an upload; the verification retries for about five minutes. `docker-image` does not
block the release and can be rerun on its own; retrying an older release never moves `MAJOR.MINOR` or
`latest` back.

Do not use **Re-run all jobs** after an upload: `qualify` then fails because the version exists. Until
the `testpypi` job uploads, a failed release can be fixed with a new commit and dispatched again with the
same version. After an upload, any change to the source or the artifacts needs a new version.

## Configuration

The workflow publishes through PyPI trusted publishing and needs no stored tokens for the indexes. The
repository settings it depends on are:

- trusted publishers for `sb-ncbr/chargefw`, workflow `release.yml`, on TestPyPI (environment
  `testpypi`) and PyPI (environment `pypi`);
- GitHub environments `testpypi`, `pypi` with the maintainer as required reviewer, and `dockerhub` with
  the `DOCKERHUB_TOKEN` secret, a Read & Write access token of the `sbncbr` Docker Hub account; all three
  limited to `main`;
- the Read the Docs automation rule that activates tags matching `^v\d+\.\d+\.\d+`;
- the Figshare GitHub integration, which archives every GitHub Release.
