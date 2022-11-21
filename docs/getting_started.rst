.. _getting started:
===============
Getting Started
===============

Public Datasets
===============

`TNBC1 data`_: a triple-negative breast cancer (TNBC) sample that was assayed by droplet-based scRNA-seq (10x Genomics), from `Gao et al. (2021)`_.

`BCH869 data`_: a glioma sample BCH869 with histone H3 lysine27-to-methionine mutations (H3K27M-glioma), where 489 malignant cells and 3 non-tumour cells were probed by smart-seq2, from `Filbin et al. (2018)`_.

Raw Data matrix (in anndata format) for RDR module and BAF moudle can be downloaded and imported directly (see notebook tutorials). Data download at `demo_data`_.

Examples of XClone and steps for reproducible results are provided in Jupyter Notebook under `examples`_ folder. 
For start, please refer to tutorials analyzing `TNBC1`_ and `BCH869`_ datasets, or `TNBC1_tutorial`_ and `BCH869_tutorial`_ records step by step.

.. _examples: https://connecthkuhk-my.sharepoint.com/:f:/g/personal/rthuang_connect_hku_hk/EhB6wYPgnL1MlUGP5sLHOhQBpLv3EFG4kToa0eY7sMZDLw?e=bypaf5
.. _TNBC1: ./TNBC1_XClone_demo_v3.html
.. _BCH869: ./BCH869_XClone_demo_v3.html
.. _TNBC1_tutorial: ./TNBC1_XClone_demo_v2.html
.. _BCH869_tutorial: ./BCH869_XClone_demo_v2.html
.. _demo_data: https://connecthkuhk-my.sharepoint.com/:f:/g/personal/rthuang_connect_hku_hk/EnKri0rS-ZpHl0VGVHUp4k0B_3iZ_gpD-obVuDwEMQUieQ?e=k0eR4T
.. _TNBC1 data: https://connecthkuhk-my.sharepoint.com/:f:/g/personal/rthuang_connect_hku_hk/Etlhi3gMu_VJuhmtrQiQRO4BRu4VVxIE_yL3Mt6iQ10kkA?e=zV0qbe
.. _BCH869 data: https://connecthkuhk-my.sharepoint.com/:f:/g/personal/rthuang_connect_hku_hk/EhnxMmkOFsNOto8XN0OYNr0BNVAvZOem3SKFcpjBKMTJFw?e=0e73Rg
.. _Gao et al. (2021): https://www.nature.com/articles/s41587-020-00795-2
.. _Filbin et al. (2018): DOI: 10.1126/science.aao4750



XClone on GX109
===============

XClone also provides integrated functions (RDR BAF and Combination) for CNV analysis by default 
whilst specific configurations might need to be adjusted accordingly.






