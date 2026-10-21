# Sphinx configuration for the nf-MiTo documentation

project = "nf-MiTo"
author = "Andrea Cossa"
copyright = "2025, Andrea Cossa"
release = "2.0.0"

extensions = ["myst_parser"]
myst_enable_extensions = ["colon_fence", "deflist"]
myst_heading_anchors = 3

source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
master_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = []
html_theme_options = {"collapse_navigation": False, "navigation_depth": 2}
