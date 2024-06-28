[![doi 10.5281/zenodo.4915837](https://zenodo.org/badge/DOI/10.5281/zenodo.4915837.svg)](https://doi.org/10.5281/zenodo.4915837)

# Optimising expression quantitative trait locus mapping workflows for single-cell studies

This repository contains scripts for the analysis of both real and simulated single-cell RNA-seq data for our paper:

Cuomo*, Alvari*, Azodi* et al, [Optimizing expression quantitative trait locus mapping workflows for single-cell studies](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02407-x), Genome Biology (2021).

## Analysis Scripts

The following folders contain scripts to reproduce the figures and analyses presented in the paper.
In particular:

* [Plotting Notebooks](ipsc) contains notebooks to reproduce our figures (both main and supplementary).

* [Simulation Scripts](simulations) contains scripts to estimate parameters from real data and then set up simulations). For an introduction to splatPop for simulating population scale single-cell RNA-seq data, see the [vignette](http://www.bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatPop.html). More detail about splatPop cen be found in our [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02546-1)!

## Software Requirements

- splatter v1.17.1
- VariantAnnotation v1.36.0
- tidyverse v1.3.1
- SingleCellExperiment v1.12.0
- fitdistrplus v1.1-3
- scater v1.18.6
- data.table v1.14.0
