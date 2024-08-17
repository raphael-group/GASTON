Installation
============


PyPI
----

You can directly install GASTON via pip::



    pip install gaston-spatial



Alternatively, you can directly install the conda environment from the `environment.yml` file. First install conda environment from `environment.yml` file::


    cd GASTON
    conda env create -f environment.yml


Then install GASTON using pip (will add to pypi soon!)::

    conda activate gaston-package
    pip install -e .

Conda
-----

Coming soon!


Development version
-------------------
To install the development version, clone the repository and install using pip::

    pip install git+https://github.com/raphael-group/gaston@main
