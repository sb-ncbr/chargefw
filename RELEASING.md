# ChargeFW Maintainer Release Procedure

A release is a reviewed version commit on `main`, validated locally, and published by the
[release workflow](.github/workflows/release.yml). The workflow builds the artifacts once, publishes
them to TestPyPI, tags the commit, and publishes the same files to PyPI after maintainer approval. It
then creates the GitHub Release and pushes the CLI image to Docker Hub.

## Current wheel matrix

The release target is one Linux x86-64 wheel built against CPython's stable ABI (`cp312-abi3`), which
serves CPython 3.12 and later. `cibuildwheel` builds it once in a manylinux container, repairs its native
libraries, and installs it and runs the Python test suite with CPython 3.12, 3.13, and 3.14. The wheel
carries `manylinux_2_27_x86_64` and `manylinux_2_28_x86_64` tags.

The wheel also installs on CPython 3.15, but it is not tested there, because the pinned optional Gemmi
0.7.4 package does not currently provide a usable CPython 3.15 wheel. Other operating systems and
architectures are not yet qualified.

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

CI tests only `gcc-debug`. Run the complete GCC, Clang, release, sanitizer, and clang-tidy matrix and build
the distributable CLI image on the version commit:

```bash
./utils/release/02-validate.sh
```

It requires `cmake`, `ctest`, `pre-commit`, the development dependencies from
[DEVELOPMENT.md](DEVELOPMENT.md), and Podman or Docker; set `CHARGEFW_CONTAINER_ENGINE=docker` to choose
Docker when both are installed. The sanitizer builds run sequentially with one build job.

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
| `build` | Builds the sdist and the wheels from it, checks the artifact set, notices, and metadata, runs a smoke calculation, and stores the artifacts with a `SHA256SUMS` manifest for 90 days |
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
