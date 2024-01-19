====================
XClone preprocessing
====================

XClone takes three cell by gene matrices as input: allele-specific AD and DP matrices and 
the total read depth matrix.
XClone preprocessing pipeline is aimed to generate the three matrices from **SAM/BAM/CRAM** files.
We recommend you use xcltk_ as the `preprocessing pipeline`_. Before that, you need :ref:`prepare the data <Preparing data>`.

* Step1: :ref:`Install xcltk <xcltk installation>` and run xctlkt :ref:`RDR module <xcltk RDR>` and :ref:`BAF module <xcltk BAF>`, independently.
* Step2: :ref:`Load xctlkt RDR module <rdr load>`.
* Step3: :ref:`Load xctlkt BAF module <baf load>`.
* Step4: XClone analysis, see page :ref:`Getting started <getting started>`.


.. _xcltk installation:

Tool installation
=================

Preprocessing via xcltk 
-----------------------

xcltk is a toolkit for XClone count generation. 
We recommend to use xcltk_ Read depth count matrix and allelic count matrix.
xcltk is avaliable through pypi. To install, type the following command line, and add -U for upgrading ::

    pip install -U xcltk

Alternatively, you can install from this GitHub repository for latest (often development) version by following command line ::

    pip install -U git+https://github.com/hxj5/xcltk


Required data
=============

.. _rdr load:

XClone RDR module
-----------------
For RDR module, we need 3 files to create ``Anndata`` format data with layer ``raw_expr``, 
cell annotation in ``obs`` and feature annotation in ``var``.

* RDR matrix
* cell annotation
* feature annotation

**load RDR demo data**

.. code-block:: python

   import xclone
   RDR_adata = xclone.data.tnbc1_rdr()
   ## preview data details
   RDR_adata
   RDR_adata.obs
   RDR_adata.var


**How to build anndata from the files**

Specify path of the files, and use ``xclone.pp.xclonedata`` to get the anndata format.
``mtx_barcodes_file`` at least include cell barcodes ID. If there are more cell annotations
in other file, can use ``xclone.pp.extra_anno`` to import.

.. code-block:: python

   import xclone
   data_dir = "xxx/xxx/xxx/"
   RDR_file = data_dir + "xcltk.rdr.mtx" 
   mtx_barcodes_file = data_dir + "barcodes.tsv" # cell barcodes
   regions_anno_file = data_dir + "features.tsv" # feature annnotation
   xclone.pp.xclonedata(RDR_file, 
                        data_mode = 'RDR', 
                        mtx_barcodes_file, 
                        regions_anno_file, 
                        genome_mode = "hg38_genes", 
                        data_notes = None)
   
   RDR_adata = xclone.pp.extra_anno(RDR_adata, anno_file, barcodes_key = "cell",
               cell_anno_key = ["Clone_ID", "Type", "cell_type"], sep = ",")
   # default sep = ",", also support "\t"

if you use the ``xcltk`` tool to prepare the input matrix, then you could find it easier to 
use default feature annotation after you specify the genome_mode (include: "hg19_genes", 
"hg38_genes" and also default 5M length blocks annotation).

.. code-block:: python

   RDR_adata = xclone.pp.xclonedata(RDR_file, 'RDR', mtx_barcodes_file, genome_mode = "hg19_genes")

.. _baf load:

XClone BAF module
-----------------

For BAF module, we need 4 files to create ``Anndata`` format data with layers ``AD`` and ``DP``, 
cell annotation in ``obs`` and feature annotation in ``var``.

* AD matrix
* DP matrix
* cell annotation
* feature annotation

**load BAF demo data**

.. code-block:: python

   import xclone
   BAF_adata = xclone.data.tnbc1_baf()
   ## preview data details
   BAF_adata
   BAF_adata.obs
   BAF_adata.var


**How to build anndata from the files**

Specify path of the files, and use ``xclone.pp.xclonedata`` to get the anndata format, similar with 
RDR module. Here the ``AD_file`` and ``DP_file`` are sparse matrix imported as ``AD`` and ``DP`` layers.

.. code-block:: python

   import xclone
   data_dir = "xxx/xxx/xxx/"
   AD_file = data_dir + "AD.mtx"
   DP_file = data_dir + "DP.mtx"
   mtx_barcodes_file = data_dir + "barcodes.tsv" # cell barcodes
   # use default gene annotation
   BAF_adata = xclone.pp.xclonedata([AD_file, DP_file], 'BAF', 
                                    mtx_barcodes_file, 
                                    genome_mode = "hg19_genes")
   BAF_adata = xclone.pp.extra_anno(BAF_adata, anno_file, barcodes_key = "cell",
               cell_anno_key = ["Clone_ID", "Type", "cell_type"], sep = ",")

.. _Preparing data:

Preparing data
==============

Detail instructions on how to prepare the data for generating Anndata for RDR module and
BAF module. 
Both part need annotation data for cell and genome features. We recommend you prepare the 
annotation data as follows.

Annotation data
---------------

**Feature annotation**

Feature annotation at least includes ``chr``, ``start``, ``stop``, ``arm`` information and 
in chr1-22,X,Y order for intuitive visualization and analysis. Here are two feature annotation
examples in `XClone` and you can load as your annotation file. If you use xcltk_ pipeline, there 
are default annotations provided.

.. code-block:: python

   import xclone
   hg38_genes = xclone.pp.load_anno(genome_mode = "hg38_genes")
   hg38_blocks = xclone.pp.load_anno(genome_mode = "hg38_blocks")


.. csv-table:: Feature (genes) annotation sample in hg38
   :file: ./tutorial_data/hg38_genes_sample.csv
   :widths: 20, 20, 10, 10, 10, 10, 10, 10
   :header-rows: 1

.. csv-table:: Feature (blocks) annotation sample in hg38
   :file: ./tutorial_data/hg38_blocks_sample.csv
   :widths: 30, 30, 20, 20
   :header-rows: 1

**Cell annotation**

* cell barcodes

`barcodes_file` include barcodes without any hearder.

.. csv-table:: barcodes_sample
   :file: ./tutorial_data/barcodes_sample.tsv
   :widths: 100
   :header-rows: 0

* cell annotation

Cell annotation (`anno_file`) at least includes ``cell``, ``cell_type``
information (Tumor or Normal, T/N), where ``cell`` is the key of cell barcodes.

.. csv-table:: cell annotation sample
   :file: ./tutorial_data/cell_anno_sample.csv
   :widths: 20, 20, 20, 10, 10, 10, 10
   :header-rows: 1


Prepare the allele-specific data (BAF) and expression data (RDR)
----------------------------------------------------------------
XClone takes 2 cell by features (genes/blocks) integer allelic AD and DP count matrices as BAF input, and it takes a cell by features (genes/blocks) integer UMI/read count matrix as RDR input. 
For BAF, we recommend using ``xcltk`` tool to get the two allelic AD and DP matrices. For RDR, you may use ``xcltk``, ``10x CellRanger`` or any other expression quantification tools to get the RDR UMI/read count matrix.

See `xcltk_preprocess`_ for details of how to prepare BAF and RDR data.


.. _xcltk: https://pypi.org/project/xcltk/
.. _preprocessing pipeline: https://github.com/hxj5/xcltk/tree/master/preprocess
.. _xcltk_preprocess: https://github.com/hxj5/xcltk/tree/master/preprocess

