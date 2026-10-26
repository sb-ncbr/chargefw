"""Sphinx configuration for the ChargeFW documentation."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "_ext"))

project = "ChargeFW"
author = "ChargeFW contributors"

extensions = ["myst_parser", "parameter_sets"]
source_suffix = {".md": "markdown"}
myst_heading_anchors = 3
myst_enable_extensions = ["dollarmath", "amsmath"]
exclude_patterns = ["_build", "_generated"]

html_theme = "furo"
html_title = "ChargeFW"
