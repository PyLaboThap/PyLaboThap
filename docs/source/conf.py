# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PyLaboThap'
copyright = '2024, Basile Chaudoire, Elise Neven'
author = 'Elise Neven, Basile Chaudoire'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc'] # https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html to generate documentation from docstrings

templates_path = ['_templates']
exclude_patterns = []

import os
import sys
# Update this path to the directory containing your library source code
sys.path.insert(0, os.path.abspath('../library'))


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


