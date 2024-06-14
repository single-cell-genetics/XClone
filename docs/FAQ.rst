FAQ
===

.. contents:: Contents
   :depth: 2
   :local:

Top Questions
-------------

Input and Output
----------------

CNV Modules
-----------

Parameters
----------

Version 0.3.4 specific issues
-----------------------------

FileNotFoundError
~~~~~~~~~~~~~~~~~


You may encounter a `FileNotFoundError` like the one shown below:

.. code-block:: python

    FileNotFoundError                         Traceback (most recent call last)
    /tmp/pbs.1280697.xomics/ipykernel_79995/464421726.py in <module>
          9     mtx_barcodes_file,
         10     genome_mode = "hg38_genes",
    ---> 11     data_notes = None
         12 )
         13 

    ~/anaconda3/envs/xclone0.3.4/lib/python3.7/site-packages/xclone/preprocessing/_data.py in xclonedata(X, data_mode, mtx_barcodes_file, regions_anno_file, genome_mode, data_notes)
        232     ### var anno
        233     if regions_anno_file is None:
    --> 234         regions_anno = load_anno(genome_mode)
        235     else:
        236         regions_anno = pd.read_table(regions_anno_file, header = None, index_col=0)

    ~/anaconda3/envs/xclone0.3.4/lib/python3.7/site-packages/xclone/preprocessing/_anno_data.py in load_anno(genome_mode)
         21     # stream.read()
         22     if genome_mode == "hg38_genes":
    --> 23         stream = pkg_resources.resource_stream(__name__, '../data/anno_data/annotate_genes_hg38_update.txt')
         24     if genome_mode == "hg38_blocks":
         25         stream = pkg_resources.resource_stream(__name__, '../data/anno_data/annotate_blocks_hg38_update.txt')

    ~/anaconda3/envs/xclone0.3.4/lib/python3.7/site-packages/pkg_resources/__init__.py in resource_stream(self, package_or_requirement, resource_name)
       1159         """Return a readable file-like object for specified resource"""
       1160         return get_provider(package_or_requirement).get_resource_stream(
    -> 1161             self, resource_name
       1162         )
       1163 

    ~/anaconda3/envs/xclone0.3.4/lib/python3.7/site-packages/pkg_resources/__init__.py in get_resource_stream(self, manager, resource_name)
       1630 
       1631     def get_resource_stream(self, manager, resource_name):
    -> 1632         return open(self._fn(self.module_path, resource_name), 'rb')
       1633 
       1634     def _get(self, path):

    FileNotFoundError: [Errno 2] No such file or directory: '/home/rthuang/anaconda3/envs/xclone0.3.4/lib/python3.7/site-packages/xclone/preprocessing/../data/anno_data/annotate_genes_hg38_update.txt'

You may download the `anno_data` from the following URL and place the files under `/data/anno_data`:

`https://github.com/single-cell-genetics/XClone/tree/master/xclone/data/anno_data`


Version 0.3.5 specific issues
-----------------------------

Not reproted so far.


Dependency
----------

ModuleNotFoundError: No module named 'requests'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may encounter an error indicating that the `requests` module is not found (in v0.3.4, v0.3.5). To resolve this, you can install the package manually:

.. code-block:: bash

    pip install requests

ModuleNotFoundError: No module named 'importlib.metadata'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may encounter a `ModuleNotFoundError` indicating that the `importlib.metadata` module is not found (in v0.3.4, v0.3.5) like the one shown below:

.. code-block:: python

   ModuleNotFoundError                       Traceback (most recent call last)
   /tmp/pbs.1280697.xomics/ipykernel_11066/2968024211.py in <module>
   ----> 1 RDR_Xdata = xclone.model.run_RDR(RDR_adata, config_file = xconfig)

   ~/anaconda3/envs/xclone0.3.5/lib/python3.7/site-packages/xclone/model/xclone_rdr_wrap.py in run_RDR(RDR_adata, verbose, run_verbose, config_file)
      225                                               low_dim=False, run_KNN=True,
      226                                               KNN_neighbors = KNN_neighbors,
   --> 227                                               copy=True)
      228 
      229     if multi_refcelltype:

   ~/anaconda3/envs/xclone0.3.5/lib/python3.7/site-packages/xclone/model/_RDR_process.py in extra_preprocess(adata, ref_celltype, cluster_key, avg_key, depth_key, low_dim, run_KNN, KNN_neighbors, copy)
      78         adata.X = np.log(adata.layers['ref_normalized'] + 0.3)
      79         sc.pp.pca(adata)
   ---> 80         sc.pp.neighbors(adata, n_neighbors = KNN_neighbors, n_pcs=40)
      81         ## Notes: connectivities and distances can be slightly different every run
      82         ## even the random_state = 0 (default).

   ~/anaconda3/envs/xclone0.3.5/lib/python3.7/site-packages/scanpy/neighbors/__init__.py in neighbors(adata, n_neighbors, n_pcs, use_rep, knn, random_state, method, metric, metric_kwds, key_added, copy)
      145         metric=metric,
      146         metric_kwds=metric_kwds,
   --> 147         random_state=random_state,
      148     )
      149 

   ~/anaconda3/envs/xclone0.3.5/lib/python3.7/site-packages/scanpy/neighbors/__init__.py in compute_neighbors(self, n_neighbors, knn, n_pcs, use_rep, method, random_state, write_knn_indices, metric, metric_kwds)
      813                 knn_distances,
      814                 self._adata.shape[0],
   --> 815                 self.n_neighbors,
      816             )
      817         # overwrite the umap connectivities if method is 'gauss'

   ~/anaconda3/envs/xclone0.3.5/lib/python3.7/site-packages/scanpy/neighbors/__init__.py in _compute_connectivities_umap(knn_indices, knn_dists, n_obs, n_neighbors, set_op_mix_ratio, local_connectivity)
      390         # umap 0.5.0
      391         warnings.filterwarnings("ignore", message=r"Tensorflow not installed")
   --> 392         from umap.umap_ import fuzzy_simplicial_set
      393 
      394     X = coo_matrix(([], ([], [])), shape=(n_obs, 1))

   ~/anaconda3/envs/xclone0.3.5/lib/python3.7/site-packages/umap/__init__.py in <module>
      34 import numba
      35 
   ---> 36 from importlib.metadata import version, PackageNotFoundError
      37 
      38 try:

   ModuleNotFoundError: No module named 'importlib.metadata'

To resolve this, you can install the package manually:

.. code-block:: bash

    pip install importlib-metadata


If the problem still exists, you can check
.. code-block:: bash

    pip show importlib-metadata

and will get the information

.. code-block:: bash

    Name: importlib-metadata
    Version: 6.7.0
    Summary: Read metadata from Python packages
    Home-page: https://github.com/python/importlib_metadata
    Author: Jason R. Coombs
    Author-email: jaraco@jaraco.com
    License: 
    Location: /home/rthuang/anaconda3/envs/xclone0.3.4/lib/python3.7/site-packages
    Requires: typing-extensions, zipp
    Required-by: anndata, numba, pynndescent, scanpy

And check if you can pip install the packages it required by again. Here we tested reinstall scanpy and numba, then it works.
The most import step you may try is:

.. code-block:: bash

    pip install scanpy


