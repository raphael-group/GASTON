# GASTON

GASTON is an interpretable deep learning model for learning a _topographic map_ of a tissue slice from spatially resolved transcriptomics (SRT) data. Specifically, GASTON models gene expression topography by learning the _isodepth_, a 1-D coordinate that smoothly varies across a tissue slice.

<p align="center">
<img src="https://github.com/raphael-group/GASTON/blob/main/docs/_static/img/gaston_figure-github.png?raw=true" height=600/>
</p>

The isodepth and topographic map learned by GASTON have many downstream applications including
* identifying **spatial domains**,
* inferring **continuous gradients** and **discontinuous jumps** in gene expression,
* deriving maps of spatial variation in **cell type organization**,
* analyzing continuous gradients in the **tumor microenvironment**

<p align="center">
  <img src="https://github.com/raphael-group/GASTON/blob/main/docs/_static/img/gaston_figure-github2.png?raw=true" height=500/>
</p>

## Installation
First install conda environment from `environment.yml` file:

```
cd GASTON
conda env create -f environment.yml
```

Then install GASTON using pip.

```
conda activate gaston-package
pip install -e .
```

We will add GASTON to PyPI and Bioconda soon!

## Getting started

Check out our [readthedocs](https://gaston.readthedocs.io/en/latest/index.html), which includes tutorials for the cerebellum and tumor analyses. GASTON requires an NxG gene expression matrix (e.g. UMI counts) and an Nx2 spatial location matrix, which can be provided or read from Space Ranger output.

The data to run the tutorials is located in `docs/notebooks/tutorials/`. Note that due to Github size constraints, you have to download the counts matrices for both analyses from [Google Drive](https://drive.google.com/drive/folders/1GiibZwhpzlur8C1hNHa1g7I4jsc1Gmn7?usp=sharing). 

Olfactory bulb/MERFISH tutorials coming soon!

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
* scanpy
* squidpy

See full list in `environment.yml` file.

## Citations

The GASTON manuscript is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.10.10.561757v1). If you use GASTON for your work, please cite our paper.

```
@article {Chitra2023,
	author = {Uthsav Chitra and Brian J. Arnold and Hirak Sarkar and Cong Ma and Sereno Lopez-Darwin and Kohei Sanno and Benjamin J. Raphael},
	title = {Mapping the topography of spatial gene expression with interpretable deep learning},
	elocation-id = {2023.10.10.561757},
	year = {2023},
	doi = {10.1101/2023.10.10.561757},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2023/10/13/2023.10.10.561757},
	eprint = {https://www.biorxiv.org/content/early/2023/10/13/2023.10.10.561757.full.pdf},
	journal = {bioRxiv}
}


```
