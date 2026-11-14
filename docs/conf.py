"""Sphinx configuration for the ChargeFW documentation."""

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "_ext"))

SCHEMAS = Path(__file__).parent.parent / "schemas"

project = "ChargeFW"
author = "ChargeFW contributors"

extensions = ["myst_parser", "parameter_sets"]
source_suffix = {".md": "markdown"}
myst_heading_anchors = 3
myst_enable_extensions = ["dollarmath", "amsmath"]
exclude_patterns = ["_build", "_generated"]

html_theme = "furo"
html_title = "ChargeFW"


def copy_schemas(app, exception):
    """Publish the JSON schemas at schemas/<name>, the address in their $id."""
    if exception is None and app.builder.format == "html":
        shutil.copytree(SCHEMAS, Path(app.outdir) / "schemas", dirs_exist_ok=True)


def setup(app):
    app.connect("build-finished", copy_schemas)
