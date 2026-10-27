"""The chargefw command installed with the package."""

import json
import shutil
import subprocess
import sys
from pathlib import Path

WATER = Path(__file__).parents[1] / "fixtures" / "synthetic" / "sdf" / "water.sdf"


def run_command(*arguments: str) -> subprocess.CompletedProcess[str]:
    # The same entry point that the installed chargefw script calls.
    launcher = "from chargefw._cli import main; main()"
    return subprocess.run(
        [sys.executable, "-c", launcher, *arguments],
        capture_output=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )


def test_command_finds_the_bundled_parameter_sets(tmp_path: Path) -> None:
    result = run_command("calculate", str(WATER), str(tmp_path))

    assert result.returncode == 0, result.stderr
    output = json.loads((tmp_path / "water.chargefw.json").read_text(encoding="utf-8"))
    assert output["status"] == "success"


def test_command_accepts_names_outside_the_system_code_page(tmp_path: Path) -> None:
    # Windows code pages such as 1252 cannot represent these names.
    directory = tmp_path / "Tomáš 分子"
    directory.mkdir()
    shutil.copy(WATER, directory / "vodička 水.sdf")

    result = run_command("calculate", str(directory / "vodička 水.sdf"), str(directory / "výstup"))

    assert result.returncode == 0, result.stderr
    output = json.loads(
        (directory / "výstup" / "vodička 水.chargefw.json").read_text(encoding="utf-8")
    )
    assert output["status"] == "success"


def test_command_passes_on_the_exit_status(tmp_path: Path) -> None:
    result = run_command("calculate", str(tmp_path / "missing.sdf"), str(tmp_path))

    assert result.returncode != 0
