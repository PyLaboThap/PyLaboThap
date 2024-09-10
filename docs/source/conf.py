import os
import sys

# Add the path to the library directory
sys.path.insert(0, os.path.abspath('../../library'))

# -- Project information -----------------------------------------------------
project = 'PyLaboThap'
copyright = '2024, Basile Chaudoire, Elise Neven'
author = 'Elise Neven, Basile Chaudoire'
release = '1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',   # To generate documentation from docstrings
    'sphinx.ext.napoleon',  # For Google/NumPy style docstrings
    'sphinx.ext.viewcode',  # Add links to highlighted source code
    'sphinx.ext.todo',      # Support for TODO notes
    'sphinx.ext.mathjax',  # For mathematical equations
]

autodoc_member_order = 'bysource'
autodoc_typehints = 'description'

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**/.ipynb_checkpoints']

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

def setup(app):
    app.add_css_file('custom.css')

html_title = "PyLaboThap Documentation"

html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}

# -- Options for mathematical equations ------------------------------------------------

# MathJax settings
mathjax_path = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'
