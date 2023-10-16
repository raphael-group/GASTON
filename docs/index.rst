
GASTON - Mapping the topography of spatial gene expression with interpretable deep learning
===========================================================================================

GASTON is an interpretable deep learning model for learning the *gene expression topography* of a tissue slice from spatially resolved transcriptomics (SRT) data. GASTON models gene expression topography by learning the *isodepth*, a 1-D coordinate describing continuous gene expression gradients and tissue geometry (i.e. spatial domains).

.. image:: https://raw.githubusercontent.com/raphael-group/GASTON/main/docs/_static/img/gaston_figure-github.png
    :alt: GASTON model architecture
    :width: 400px
    :align: center
    

GASTON's key applications
-------------------------
- Learns 1-d coordinate that varies smoothly across tissue slice, providing *topographic map* of gene expression in the tissue slice.
- Modeling *continuous gradients* of gene expression for individual genes, e.g. gradients of metabolism in cancer
- Identifying *tissue geometry*, i.e. arrangement of spatial domains

.. image:: https://raw.githubusercontent.com/raphael-group/GASTON/main/docs/_static/img/gaston_figure-github2.png
    :alt: GASTON model architecture
    :width: 400px
    :align: center

Manuscript
----------
Please see our manuscript :cite:`GASTON` for more details.

Getting started with GASTON
---------------------------
- Browse :doc:`notebooks/tutorials/index` for a quick start guide to GASTON.
- Discuss usage and issues on `github`_.


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

.. _github: https://github.com/raphael-group/GASTON