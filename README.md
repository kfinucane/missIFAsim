
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Simulation code for thesis: Apposite Bayesian latent factormodels for imputation in metabolomics data.

This repository contains the code used to create simulated datasets for
rIFA, TGIFA, and VI-TGIFA simulation studies, and to run the simulation
studies.

## Generating simulated data

The `data_gen.R` file generates simulated data. At the bottom of the
file, specify the desired characteristics of the data, and run the whole
file. In this repository, a simulated example of the urinary
metabolomics dataset used in the true generating process has been
provided. The real dataset is not available publicly, however, this
should allow for demonstration of how the data generation process works.
The simulated datasets used in the simulation studies are provided in
`data/p_1391_datasets` for reproducibility.

## Simulation study code

The scripts used to run the simulation studies for rIFA, TGIFA, and
VI-TGIFA are provided.

### rIFA simulation study

The `rIFA_running.R` script runs the rIFA imputation, while
`rIFA_model.R` contains the rIFA model, and `rIFA_utils.R` contains
utility functions. The `IFA_running.R` script runs the IFA imputation,
while `IFA_model.R` contains the IFA model used.

### TGIFA simulation study

The `TGIFA_running.R` script runs the TGIFA imputation, while
`tGIFA_model.R` contains the tGIFA model, and `tGIFA_utils.R` contains
utility functions.

### VI-TGIFA simulation study

The `VI-TGIFA_running.R` script runs the VI-TGIFA imputation, while
`VI-TGIFA_model.R` contains the VI-TGIFA model, and `VI-TGIFA_utils.R`
contains utility functions.

Note that in order to run the VI-TGIFA code, the VIMSFA R package is
required. This can be installed as follows:

``` r
# install remotes package if not already installed
# install.packages("remotes")
# remotes::install_github("blhansen/VI-MSFA")
```
