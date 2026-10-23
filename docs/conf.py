"""Sphinx configuration for the ChargeFW documentation."""

project = "ChargeFW"
author = "ChargeFW contributors"

extensions = ["myst_parser"]
source_suffix = {".md": "markdown"}
myst_heading_anchors = 3
exclude_patterns = ["_build"]

html_theme = "furo"
html_title = "ChargeFW"
