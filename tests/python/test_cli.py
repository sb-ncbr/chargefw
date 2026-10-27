"""The chargefw command installed with the package."""

import json
import subprocess
import sys
from pathlib import Path

WATER = Path(__file__).parents[1] / "fixtures" / "synthetic" / "sdf" / "water.sdf"


def run_command(*arguments: str) -> subprocess.CompletedProcess[str]:
    # The same entry point that the installed chargefw script calls.
    launcher = "from chargefw._cli import main; main()"
    return subprocess.run(
        [sys.executable, "-c", launcher, *arguments], capture_output=True, text=True, check=False
    )


def test_command_finds_the_bundled_parameter_sets(tmp_path: Path) -> None:
    result = run_command("calculate", str(WATER), str(tmp_path))

    assert result.returncode == 0, result.stderr
    output = json.loads((tmp_path / "water.chargefw.json").read_text())
    assert output["status"] == "success"


def test_command_passes_on_the_exit_status(tmp_path: Path) -> None:
    result = run_command("calculate", str(tmp_path / "missing.sdf"), str(tmp_path))

    assert result.returncode != 0
