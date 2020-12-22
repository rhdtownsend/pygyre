# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import re

sys.path.insert(0, os.path.abspath('../../source'))

import sphinx_rtd_theme


# -- Project information -----------------------------------------------------

project = 'PyGYRE'
author = 'Rich Townsend & The PyGYRE Team'
version = "1.1.5"
release = "1.1.5"
branch = "main"
copyright = '2020, Rich Townsend & The PyGYRE Team'


# -- General configuration ---------------------------------------------------

# Numbered figures
numfig = True

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.mathjax',
    'sphinx.ext.extlinks',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx-prompt',
    'sphinx_substitution_extensions',
]


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']


# -- Additional configuration ------------------------------------------------

# sphinx_rtd options
html_theme_options = {
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
    'logo_only': True
}

# Set master doc
master_doc = 'index'

# Set logo
html_logo = 'pygyre-logo.png'

# Set up Extlinks
extlinks = {
    'wiki': ('https://en.wikipedia.org/wiki/%s', ''),
    'ads': ('https://ui.adsabs.harvard.edu/abs/%s/abstract', ''),
    'netlib': ('https://www.netlib.org/%s', ''),
    'git': ('https://github.com/%s', ''),
    'repo': ('https://github.com/rhdtownsend/pygyre/blob/{:s}/%s'.format(branch), '')
}

# Set site-wide targets
targets = {
    'github-tarball': 'https://codeload.github.com/rhdtownsend/pygyre/tar.gz/v{0:s}'.format(release),
}

rst_prolog = '\n'.join(['.. _{:s}: {:s}'.format(x, targets[x]) for x in targets])

# Add substitutions for sphinx_substitution_extensions

rep_exts = {"release": release,
            "author": author}

for rep_ext_key, rep_ext_val in rep_exts.items():
    rst_prolog += "\n.. |{:s}| replace:: {:s}".format(rep_ext_key, rep_ext_val)

# Enable email obfuscation
email_automode = True

# Intersphinx mappig
intersphinx_mapping = {
   'gyre': ('https://gyre.readthedocs.io/en/latest/', None),
   'astropy': ('https://docs.astropy.org/en/stable/', None)
   }

# Remove module docstring from autodoc
def remove_module_docstring(app, what, name, obj, options, lines):
    if what == "module" and name == "pygyre":
        del lines[:]

def setup(app):
    app.connect("autodoc-process-docstring", remove_module_docstring)
    
