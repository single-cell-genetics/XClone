.. _api:

=================
API Documentation
=================


.. contents:: Contents
   :depth: 3
   :local:


.. automodule:: xclone
   :members:



Import XClone as::

   import xclone


Commands
========

* ``XCloneRDR``: see manual_
* ``XCloneBAF``: see manual_
* ``XCloneCombine``: see manual_

.. _manual: https://xclone-cnv.readthedocs.io/en/latest//manual.html



Functions
=========


Read / Load
-----------

.. autofunction:: xclone.preprocessing.xclonedata


.. autofunction:: xclone.preprocessing.readrdr_mtx
    

.. autofunction:: xclone.preprocessing.extra_anno



Config module-specific parmas
-----------------------------

Main config manager
~~~~~~~~~~~~~~~~~~~

.. autoclass:: xclone.PreprocessingConfig
    :members:
    :undoc-members:
    :show-inheritance:


.. autoclass:: xclone.XCloneConfig
    :members:
    :undoc-members:
    :show-inheritance:

Module specific config manager 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: xclone.XCloneGeneral_config
    :members:
    :undoc-members:
    :show-inheritance:


.. autoclass:: xclone.RDR_General_config
    :members:
    :undoc-members:
    :show-inheritance:


.. autoclass:: xclone.BAF_General_config
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: xclone.Combine_General_config
    :members:
    :undoc-members:
    :show-inheritance:


.. autoclass:: xclone.HMM_Configs
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: xclone.Smartseq_Config
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: xclone.Spatial_Config
    :members:
    :undoc-members:
    :show-inheritance:



RDR module
----------

.. autofunction:: xclone.model.run_RDR



BAF module
----------

.. autofunction:: xclone.model.run_BAF



Combine module
--------------

.. autofunction:: xclone.model.run_combine




Plotting
--------

Basic Heatmap plot
~~~~~~~~~~~~~~~~~~

.. autofunction:: xclone.plot.CNV_visualization

.. autofunction:: xclone.plot.BAF_CNV_visualization

.. autofunction:: xclone.plot.Combine_CNV_visualization




Complex Annotated Heatmap plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: xclone.plot.Complex_CNV_visualization

.. autofunction:: xclone.plot.Complex_BAF_CNV_visualization

.. autofunction:: xclone.plot.Complex_Combine_CNV_visualization





Allele frequency plot
~~~~~~~~~~~~~~~~~~~~~


.. autofunction:: xclone.plot.visualize_cell_BAF



Read Depth Ratio plot
~~~~~~~~~~~~~~~~~~~~~






VCF processing
--------------





XClone Object
-------------


Clustering
~~~~~~~~~~



Spatial analysis
~~~~~~~~~~~~~~~~



Extract Output
--------------

