import os, sys

# Import version following relative path ../ksrates/ksrates/_version.py
p = os.path.abspath('..')
sys.path.insert(1, p)
from ksrates._version import __version__

project = "ksrates"
author = "Cecilia Sensalari, Steven Maere and Rolf Lohaus"

version = __version__
release = __version__

html_theme = "sphinx_rtd_theme"

pygments_style = "sphinx"

source_suffix = '.rst'
master_doc = 'index'

# To use light grey background and no frame on code blocks
# To avoid extra white pages to force starting on the left page (like it was a book)
# To use inconsolata monospace font (looks better than Courier)
latex_elements = {

    'sphinxsetup' : 'VerbatimColor={rgb}{0.97,0.97,0.97}, VerbatimBorderColor={rgb}{1,1,1}',
    'extraclassoptions': 'openany,oneside',
    'preamble': r'''
       \usepackage{inconsolata}
    '''
}

# Add kstates white logo
html_static_path = ['_static']
html_css_files = ["custom.css"]
html_logo = "_static/logo_ksrates_long_inverted.svg"
html_theme_options = {
    'logo_only': True # Show only logo and not package name
}

# Add copy-button in code blocks
extensions = [
    'sphinx_copybutton'
]

# For local export to PDF:
# import rst2pdf
# extensions = [
#     'rst2pdf.pdfbuilder',
# ]

# COMMAND TO BUILD THE DOCUMENTATION LOCALLY
# Install v6.2.1. due to bug in >= 7.0.0 with "sphinx_rtd_theme" theme [to date: March2024]
# pip install -U sphinx==6.2.1
# For LaTeX: sphinx-build -M latexpdf source/ abs/path/outdir
# It will generate the PDF version in outdir/latex
# For HTML: sphinx-build -b html docs/ docs/outdir