# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
import time
sys.path.insert(0, os.path.abspath('../../femmt'))


# -- Project information -----------------------------------------------------

project = 'FEM Magnetics Toolbox'
copyright = '{}, {}'.format(time.strftime('%Y'),  'LEA, Paderborn University')
author = 'LEA-UPB'

# The full version, including alpha/beta/rc tags
release = '0.1.2'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
sys.path.insert(0, os.path.abspath('_extensions'))
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx.ext.todo', 'sphinx_multiversion']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates/']

# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = None
todo_include_todos = True

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
smv_released_pattern = r'^tags/.*$'

html_theme_options = {
    'logo_only': False,
    'display_version': True,
    'collapse_navigation': False,
}

# replace "view page source" with "edit on github" in Read The Docs theme
#  * https://github.com/readthedocs/sphinx_rtd_theme/issues/529

html_context = {
    'display_github': True,
    'github_user': 'upb-lea',
    'github_repo': 'FEM_Magnetics_Toolbox',
    'github_version': 'main/sphinx/source/',
}

# Code for adding html side bars to theme
html_sidebars = {
    '**': [
        'versioning.html',
    ],
 }
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']
