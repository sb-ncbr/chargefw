"""The chargefw command: starts the command line installed beside the extension."""

import os
import subprocess
import sys
from pathlib import Path
from typing import NoReturn


def main() -> NoReturn:
    name = "_chargefw_cli.exe" if sys.platform == "win32" else "_chargefw_cli"
    executable = str(Path(__file__).with_name(name))
    arguments = ["chargefw", *sys.argv[1:]]
    if sys.platform == "win32":
        # Windows cannot replace the running process; wait for the program and pass on its status.
        raise SystemExit(subprocess.call(arguments, executable=executable))
    os.execv(executable, arguments)
