

GASTON - Mapping the topography of spatial gene expression with interpretable deep learning
===========================================================================================

Overview
--------

GASTON is an interpretable deep learning model for learning the _topography_ of a tissue slice, using spatially resolved transcriptomics (SRT) data. Specifically, GASTON models gene expression topography by learning the _isodepth_, a 1-D coordinate describing tissue geometry (i.e. spatial domains) and continuous gene expression gradients.

Installation
------------

First install conda environment from `environment.yml` file:


    cd GASTON
    conda env create -f environment.yml


Then install GASTON using pip (will add to pypi soon!)

    conda activate gaston-package
    pip install -e .



Software dependencies
---------------------
- torch
- matplotlib
- pandas
- scikit-learn
- numpy
- jupyterlab
- seaborn
- tqdm
- scipy
- scanpy

See the `environment.yml` file

## Getting started
Try out the Jupyter notebook tutorial: `tutorial.ipynb`. Note you have to download the counts matrix here: https://drive.google.com/drive/folders/1GiibZwhpzlur8C1hNHa1g7I4jsc1Gmn7?usp=sharing

A readthedocs is coming soon!
