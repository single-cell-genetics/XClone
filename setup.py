"""
xclone - Inference of clonal Copy Number Variation in single cells
See: https://github.com/single-cell-genetics/XClone
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Set __version__ for the project.
exec(open("./xclone/version.py").read())

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

reqs = ['numpy>=1.18', 'scipy', 'matplotlib', 'anndata',
        'statsmodels', 'scanpy', 'h5py', 'palettable']

# seaborn

setup(
    name='xclone',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,

    description='Inference of clonal Copy Number Variation in single cells',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/single-cell-genetics/XClone',

    # Author details
    author='Rongting Huang',
    author_email='rthuang@connect.hku.hk',

    # Choose your license
    license='Apache-2.0',

    # What does your project relate to?
    keywords=['Copy Number Variation', 'allelic expression', 
              'single-cell RNA-seq', 'variantional inference'],
    
    # package data
    include_package_data=True,
    package_data={'xclone': ['data/anno_data/*.txt']},

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),

    # entry_points={
    #       'console_scripts': [
    #           ],
    #       }, 

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    
    install_requires=reqs,

    extras_require={
        'docs': [
            'sphinx_bootstrap_theme']},

    py_modules = ['xclone']

    # buid the distribution: python setup.py sdist
    # upload to pypi: twine upload dist/...
)
