|PyPI| |Docs| 

.. |PyPI| image:: https://img.shields.io/pypi/v/xclone.svg
    :target: https://pypi.org/project/xclone
.. |Docs| image:: https://readthedocs.org/projects/xclone1/badge/?version=latest
   :target: https://xclone1.readthedocs.io


.. :Author: Rongting Huang
.. :Version: 0.3.0
.. :Last viewed: Oct 20, 2022

===========================================================================================================
XClone: detection of allele-specific subclonal copy number variations from single-cell transcriptomic data
===========================================================================================================

.. image:: ./image/XClone_overview_150dpi.png
   :width: 800px
   :align: center


XClone's key features
======================
* XClone has two modules of information: the read depth ratio (RDR) module and the B-allele frequency (BAF) module, where each of the modules has its own CNV states and noise models for likelihood function. 

* XClone introduces two orthogonal methods for smoothing the probabilities of assigning CNV states for each feature in each cell.

* XClone introduces a post-step for combining the CNV states computed from the RDR and BAF modules.


Latest news
============


References
==========

For details of the method, please checkout our paper :cite:`huang2023xclone`.


Support
=======

Found a bug or would like to see a feature implemented? Feel free to submit an
`issue <https://github.com/single-cell-genetics/XClone/issues/new>`_.
Have a question or would like to start a new discussion? Head over to
`GitHub discussions <https://github.com/single-cell-genetics/XClone/discussions>`_.
In either case, you can also always send us an `email <mailto:rthuang@connect.hku.hk>`_.
Your help to improve XClone is highly appreciated.
For further information visit `xclone.org <https://xclone-cnv.readthedocs.io/en/latest/>`_.


.. toctree::
   :caption: Main
   :maxdepth: 1
   :hidden:

   about
   installation
   FAQ
   release
   references

.. toctree::
   :caption: Tutorials
   :maxdepth: 1
   :hidden:
   
   tutorials
   getting_started
   preprocessing
   TNBC1_XClone_tutorials
   BCH869_XClone_tutorials
