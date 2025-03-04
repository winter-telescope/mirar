# pylint: disable=all
from mirar.paths import PACKAGE_NAME, __version__

# Configuration file for the Sphinx documentation builder.

# -- Project information

project = PACKAGE_NAME
copyright = "2022, Robert Stein"
author = "Robert Stein"

release = __version__
version = __version__

# -- General configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx_mdinclude",
]

doctest_global_setup = """
import contextlib

with contextlib.redirect_stdout(None):
    from astroquery.gaia import Gaia
"""


intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]

# -- Options for HTML output

html_theme = "sphinx_rtd_theme"

# -- Options for EPUB output
epub_show_urls = "footnote"
