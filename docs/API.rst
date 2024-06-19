.. _api:

===
API
===

.. automodule:: xclone

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

.. autofunction:: xclone.preprocessing.extra_anno



Config module-specific parmas
-----------------------------

.. autofunction:: xclone.XCloneConfig



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

