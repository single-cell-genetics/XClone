FAQ
===

.. contents:: Contents
   :depth: 2
   :local:

Top Questions
-------------

Input and Output
----------------

CNA Modules
-----------

Parameters
----------


Python Environment
------------------

For detailed information about Python environments, refer to the `Python environments issue <https://github.com/single-cell-genetics/XClone/issues/6>`_.

XClone versions 0.3.4 and 0.3.5 perform optimally on Python 3.7 (not >3.7 due to dependency issues). We recommend using conda to manage the environment. 
You can use the following environment files to create a suitable conda environment for running XClone:

- `xclone0.3.4Python3.7_env.yml <https://github.com/Rongtingting/xclone-data/blob/main/XClone_env/xclone0.3.4Python3.7_env.yml>`__
- `xclone0.3.5Python3.7_env.yml <https://github.com/Rongtingting/xclone-data/blob/main/XClone_env/xclone0.3.5Python3.7_env.yml>`__


XClone version 0.3.6 works well with Python 3.7 and Python >=3.9 (excluding 3.8 due to dependency issues).

We strongly recommend using Python >=3.9 for XClone versions >=0.3.6.

You can refer to the `xclone0.3.6Python3.9_env.yml <https://github.com/Rongtingting/xclone-data/blob/main/XClone_env/xclone0.3.6Python3.9_env.yml>`__.

Creating a Conda environment from a `.yml` file is straightforward. Below are the steps to create an environment using the provided `xclone0.3.6Python3.9_env.yml` file.

Steps to Create a Conda Environment from a `.yml` File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Download the `.yml` file:**

   First, download the `xclone0.3.6Python3.9_env.yml` file from the GitHub repository or ensure it's available locally.

2. **Open a terminal or command prompt:**

   Ensure you have Conda installed. If not, you can download it from the `Anaconda website <https://www.anaconda.com/products/distribution>`__ or use `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__ for a lighter version.

3. **Navigate to the directory containing the `.yml` file:**

   .. code-block:: sh

      cd path/to/your/yml/file

4. **Create the Conda environment:**

   Use the `conda env create` command to create an environment from the `.yml` file:

   .. code-block:: sh

      conda env create -f xclone0.3.6Python3.9_env.yml

5. **Activate the environment:**

   After the environment is created, activate it using:

   .. code-block:: sh

      conda activate xclone_advancepy39


   The name of the environment (`xclone_advancepy39`) is derived from the `name` field in the `.yml` file.

6. **Verification**
   
   After creating and activating the environment, you can verify the installation by checking the installed packages:

   .. code-block:: sh

      conda list

   This command will display all the packages and their versions installed in your Conda environment.

By following these steps, you should be able to create and activate a Conda environment based on the specifications provided in the `xclone0.3.6Python3.9_env.yml` file.


Dependency requirements
~~~~~~~~~~~~~~~~~~~~~~~~

For dependency requirements recommended by Poetry, see `xclone0.3.6Python3.9Project.toml <https://github.com/Rongtingting/xclone-data/blob/main/XClone_env/xclone0.3.6Python3.9Project.toml>`__.



XClone Version 0.3.5 specific issues
------------------------------------

Only support python == 3.7.
Recommend XClone version 0.3.6 for Python >=3.9 environment.
No other issues reproted so far.



XClone Version 0.3.4 specific issues
------------------------------------
Only support python == 3.7.
Recommend XClone version 0.3.6 for Python >=3.7 environment.

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




Dependency
----------

ModuleNotFoundError: No module named 'requests'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may encounter an error indicating that the `requests` module is not found (in v0.3.4, v0.3.5). To resolve this, you can install the package manually:

.. code-block:: bash

    pip install requests

This Dependency issues solved in XClone version >=0.3.6.


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

This Dependency issues solved in XClone version >=0.3.6 (for Python >=3.9).

For XClone version 0.3.6 (Python ==3.7), the ModuleNotFoundError: No module named 'importlib.metadata' error indicates that the importlib.metadata module is not found, 
which is unexpected given that importlib-metadata is included in setup.py and installed. This issue is likely due to the importlib.metadata module being available only in Python 3.8 and later. 
Since you are using Python 3.7, you need to install the backport package importlib-metadata.
