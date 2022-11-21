====================
XClone preprocessing
====================

XClone takes three cell by gene matrices as input: allele-specific AD and DP matrices and 
the total read depth matrix.
XClone preprocessing pipeline is aimed to generate the three matrices from **SAM/BAM/CRAM** files.
We recommend you use xcltk_ as the preprocessing pipeline. Before that, you need :ref:`prepare the data <Preparing data>`.

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

.. _xcltk RDR:

Prepare the expression data (RDR)
---------------------------------
XClone takes a cell by features (genes/blocks) integer UMI count matrix as input. 
We recommend you use ``xcltk`` tool to get the RDR UMI count matrix or directly use counts from 10x CellRanger.

RDR matrix could be generated by ``xcltk basefc`` command with proper settings.

.. _xcltk BAF:

Prepare the allele-specific data (BAF)
--------------------------------------

XClone takes 2 cell by features (genes/blocks) integer allelic AD and DP count matrices as input. 
We recommend you use ``xcltk`` tool to get the two allelic AD and DP matrices (BAF matrices) as following 
5 main steps(xcltk_details_).

* BAF Pre-Imputation
* BAF Imputation (or Phasing only)
* BAF pre-pileup
* BAF Pileup
* BAF Phasing SNP


BAF Pre-Imputation
^^^^^^^^^^^^^^^^^^

The first step is pre-imputation. Germline SNPs would be called from 
SAM/BAM/CRAM file and saved as VCF file, which would be processed and only 
heterozygous SNPs would be used as input of the next step 
`Sanger Imputation <Sanger Server_>`_.

Quite a few softwares and some auxiliary datasets should be installed or downloaded 
before running this step. The required resources are specified in a 
`configure file <baf_pre_impute config_>`_ [baf_pre_impute config]. It's recommended to make a copy of the configure 
file and then modify the new copy instead of modifying the original file. 

To run this step, use the script ``baf_pre_impute.sh`` in ``preprocess/baf_pre_impute``
directory with proper settings. All possible settings could be found with the 
``-h`` option,

.. code-block:: shell

   cd preprocess/baf_pre_impute
   ./baf_pre_impute.sh -h


BAF Imputation (or Phasing only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The VCF generated from previous step would be used for phasing or imputation on 
`Sanger Imputation Server <Sanger Server_>`_. There's a comprehensive introduction to
this step at `here <Sanger Wiki_>`_.

BAF pre-pileup
^^^^^^^^^^^^^^

The VCF file outputted by Sanger Imputation Server needs to be processed before used for pileuping. 

Similar to the Pre-Imputation step, this step also has a `configure file <baf_pre_pileup config_>`_
that specifies the required softwares and datasets. Still, Copy-and-Modify is recommended.

To run this step, use the script ``baf_pre_pileup.sh`` in ``preprocess/baf_post_impute``
directory with proper settings. All possible settings could be found with the ``-h`` option,

.. code-block:: shell

   cd preprocess/baf_post_impute
   ./baf_pre_pileup.sh -h


BAF Pileup
^^^^^^^^^^

Pileuping is aimed to extract variant information, mainly read depth of REF and 
ALT alleles, from ``SAM/BAM/CRAM`` file for each cell.

To run this step, the `configure file <baf_pileup config_>`_ should be provided and then use 
the script ``baf_pileup.sh`` in ``preprocess/baf_post_impute`` directory with proper settings. 
The script simply wraps cellsnp-lite and xcltk pileup. All possible settings could be 
found with the ``-h`` option,

.. code-block:: shell

   cd preprocess/baf_post_impute
   ./baf_pileup.sh -h

Alternatively, you could use other pileup tools for this step.

BAF Phasing SNP
^^^^^^^^^^^^^^^

The allele-specific read depth of each variant for each cell, i.e. the AD and DP 
matrices created in previous steps, should be phased into each self-defined blocks. 
The blocks could be segments of even size that distribute uniformly along the 
genome or other feature blocks.

To run this step, simply use the ``xcltk phase_snp`` command with proper settings.


.. _xcltk: https://pypi.org/project/xcltk/
.. _xcltk_details: https://github.com/hxj5/xcltk/tree/master/preprocess
.. _Sanger Server: https://imputation.sanger.ac.uk/
.. _Sanger Wiki: https://imputation.sanger.ac.uk/?instructions=1
.. _baf_pre_impute config: https://github.com/hxj5/xcltk/blob/master/preprocess/baf_pre_impute/baf_pre_impute.cfg
.. _baf_pre_pileup config: https://github.com/hxj5/xcltk/blob/master/preprocess/baf_post_impute/baf_pre_pileup.cfg
.. _baf_pileup config: https://github.com/hxj5/xcltk/blob/master/preprocess/baf_post_impute/baf_pileup.cfg