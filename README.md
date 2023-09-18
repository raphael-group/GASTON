# GASTON

## Overview

GASTON is an interpretable deep learning model for learning the _topography_ of a tissue slice, using spatially resolved transcriptomics (SRT) data. Specifically, GASTON models gene expression topography by learning the _isodepth_, a 1-D coordinate describing tissue geometry (i.e. spatial domains) and continuous gene expression gradients.

## Installation
Install using pip (will add to pypi soon)

```
cd GASTON
pip install -e .

```

## Software dependencies
* torch
* matplotlib
* pandas
* scikit-learn
* numpy
* jupyterlab
* seaborn
* tqdm
* scipy

See the `environment.yml` file

## Getting started
Try out the Jupyter notebook tutorial: `tutorial.ipynb`
A readthedocs is coming soon!
