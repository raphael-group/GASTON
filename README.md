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
GASTON is on pypi: https://pypi.org/project/gaston-spatial/.

```pip install gaston-spatial```

You can also directly install the conda environment from the `environment.yml` file:

```
cd GASTON
conda env create -f environment.yml
```

Then install GASTON using pip.

```
conda activate gaston-package
pip install -e .
```

Installation should take <10 minutes.

We will add GASTON to bioconda soon!

## Getting started

Check out our [readthedocs](https://gaston.readthedocs.io/en/latest/index.html), which includes tutorials for four analyses: 
* mouse cerebellum (Slide-SeqV2)
* human colorectal tumor (10x Visium)
* human breast cancer tumor (10x Xenium)
* mouse cortex (MERFISH).

Training GASTON by following the tutorials should take <1 hour.

GASTON requires an NxG gene expression matrix (e.g. UMI counts) and an Nx2 spatial location matrix, which can be provided or read from Space Ranger output.

The data to run the tutorials is located in `docs/notebooks/tutorials/`. Note that due to Github size constraints, you have to download the counts matrices for both analyses from [Google Drive](https://drive.google.com/drive/folders/1GiibZwhpzlur8C1hNHa1g7I4jsc1Gmn7?usp=sharing). 

## Software dependencies
* torch (=2.0.0)
* matplotlib (=3.8.0)
* pandas (=2.1.1)
* scikit-learn (=1.3.1)
* numpy (=1.23.4)
* jupyterlab (=4.0.6)
* seaborn (=0.12.2)
* tqdm (=4.66.1)
* scipy (=1.11.2)
* scanpy (=1.9.5)
* squidpy (=1.3.1)

See full list in `environment.yml` file. GASTON can be run on CPU or GPU.

## Citations

The GASTON manuscript is published at [Nature Methods](https://www.nature.com/articles/s41592-024-02503-3). If you use GASTON for your work, please cite our paper.

```

@article{Chitra2025,
	date = {2025/01/23},
	date-added = {2025-01-27 12:52:11 -0500},
	date-modified = {2025-01-27 12:52:11 -0500},
	doi = {10.1038/s41592-024-02503-3},
	id = {Chitra2025},
	isbn = {1548-7105},
	journal = {Nature Methods},
	title = {Mapping the topography of spatial gene expression with interpretable deep learning},
	url = {https://doi.org/10.1038/s41592-024-02503-3},
	year = {2025},
	bdsk-url-1 = {https://doi.org/10.1038/s41592-024-02503-3}}

```

## Support
For questions or comments, please file a Github issue and/or email Uthsav Chitra (uchitra@princeton.edu).
