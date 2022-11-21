============
Installation
============

Main Module XClone
==================

XClone requires Python 3.7 or later. 
We recommend to use Anaconda_ environment for version control and to avoid potential conflicts ::

    conda create -n xclone python=3.7
    conda activate xclone

XClone package can be conveniently installed via PyPI (for stable version) ::

    pip install xclone

or directly from GitHub repository (for development version)::

    pip install git+https://github.com/single-cell-genetics/XClone


Preprocessing tool xcltk
=========================

For recommeded preprocessing tool `xcltk`, refer to preprocessing page, :ref:`Tool installation <xcltk installation>`.

Dependencies
=============

Most required dependencies are automatically installed, e.g.

- `scanpy <https://scanpy-tutorials.readthedocs.io/>`_ for a few pre- and post-processing analysis
- `statsmodels <https://www.statsmodels.org/stable/index.html>`_ for GLM fitting
- `jupyter <https://jupyter.org/>`_ for running XClone CNV detection within notebooks

If you run into any issues or errors are raised during the installation process, feel free to contact us at GitHub_.

.. _Anaconda: https://www.anaconda.com/
.. _xcltk: https://pypi.org/project/xcltk/
.. _GitHub: https://github.com/single-cell-genetics/XClone
.. _`Getting Started`: getting_started
.. _`Prepare data and preprocessing`: preprocessing