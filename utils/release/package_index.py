"""Compare a release checksum manifest with the files a package index holds for that release."""

import argparse
import hashlib
import json
import platform
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

USER_AGENT = "chargefw-release-workflow"
INDEXES = {
    "testpypi": "https://test.pypi.org/pypi/chargefw/{version}/json",
    "pypi": "https://pypi.org/pypi/chargefw/{version}/json",
}


def fetch(url, timeout=30):
    with urlopen(Request(url, headers={"User-Agent": USER_AGENT}), timeout=timeout) as response:
        return response.read()


def remote_files(index, version):
    try:
        payload = json.loads(fetch(INDEXES[index].format(version=version)))
    except HTTPError as error:
        if error.code != 404:
            raise
        return {}
    return {entry["filename"]: entry for entry in payload["urls"]}


def read_manifest(path):
    expected = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        digest, artifact_path = line.split(maxsplit=1)
        expected[Path(artifact_path.lstrip("*")).name] = digest
    return expected


def compare(index, remote, manifest):
    """Return the manifest files missing from the index; fail on any conflicting file."""
    expected = read_manifest(manifest)
    conflicts = sorted(
        name for name, entry in remote.items() if expected.get(name) != entry["digests"]["sha256"]
    )
    if conflicts:
        raise SystemExit(
            f"{index} holds files that differ from the manifest: {', '.join(conflicts)}"
        )
    return sorted(expected.keys() - remote.keys())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["preflight", "download-wheel"])
    parser.add_argument("index", choices=sorted(INDEXES))
    parser.add_argument("version")
    parser.add_argument("manifest")
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    args = parser.parse_args()

    # Use one response throughout: consecutive requests can reach caches holding different lists.
    remote = remote_files(args.index, args.version)
    missing = compare(args.index, remote, args.manifest)
    if args.command == "preflight":
        print(f"{args.index}: {len(missing)} files to upload, no conflicts")
    elif missing:
        raise SystemExit(f"{args.index} does not yet list: {', '.join(missing)}")
    else:
        # The index hashes were compared with the manifest above; check the downloaded bytes too.
        # One stable-ABI wheel per platform serves every supported CPython; take this machine's.
        suffix = f"_{platform.machine()}.whl"
        wheels = [
            entry
            for name, entry in remote.items()
            if "-abi3-manylinux" in name and name.endswith(suffix)
        ]
        if len(wheels) != 1:
            raise SystemExit(
                f"expected one stable-ABI Linux {platform.machine()} wheel on {args.index}, "
                f"found {len(wheels)}"
            )
        content = fetch(wheels[0]["url"], timeout=120)
        if hashlib.sha256(content).hexdigest() != wheels[0]["digests"]["sha256"]:
            raise SystemExit("the downloaded wheel does not match its index SHA-256")
        path = args.output_dir / wheels[0]["filename"]
        path.write_bytes(content)
        print(path)


if __name__ == "__main__":
    main()
