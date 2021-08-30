project = "ksrates"
author = "Cecilia Sensalari, Steven Maere and Rolf Lohaus"

version = "1.1"
release = "1.1"

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

# For local export to PDF:
# import rst2pdf
# extensions = [
#     'rst2pdf.pdfbuilder',
# ]

# COMMAND TO BUILD THE DOCUMENTATION LOCALLY
# For LaTeX: sphinx-build -M latexpdf source/ abs/path/outdir
# It will generate the PDF version in outdir/latex
# For HTML: sphinx-build -b html source/ abs/path/outdir