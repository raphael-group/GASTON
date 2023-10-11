
GASTON - Mapping the topography of spatial gene expression with interpretable deep learning
===========================================================================================

GASTON is an interpretable deep learning model for learning the *topography* of a tissue slice, using spatially resolved transcriptomics (SRT) data. Specifically, GASTON models gene expression topography by learning the *isodepth*, a 1-D coordinate describing tissue geometry (i.e. spatial domains) and continuous gene expression gradients.

.. image:: https://raw.githubusercontent.com/raphael-group/GASTON/main/docs/_static/img/gaston_figure-V7.pdf
    
    :width: 400px
    :align: center
    :alt: GASTON model architecture

Manuscript
----------
Please see our manuscript :cite:`GASTON` for more details.

GASTON's key applications
-------------------------
- Learns 1-d coordinates that are highly correlated with the spatial coordinates of the tissue slice, and can be used to visualize the topography of gene expression in the tissue slice.
- Able to continuous gradient of gene expression for each gene, which can be used to identify genes with spatially varying expression patterns.
- Can be used to identify genes with spatially varying expression patterns, and to identify genes that are spatially co-expressed.
- Identify spatial domains in the tissue slice akin to celltypes wihch has continuous gene expression gradients.
- Plot topographical contours for *any* spatial transcriptomics dataset. 

Getting started with GASTON
---------------------------
- Browse :doc:`notebooks/tutorials/index` for a quick start guide to GASTON.
- Discuss usage and issues on `github `_.


.. toctree::
    :caption: General
    :maxdepth: 2
    :hidden:

    installation
    api
    classes
    references

.. toctree::
    :caption: Tutorial
    :maxdepth: 2
    :hidden:

    notebooks/tutorials/index