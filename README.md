# Penalized Regression Calibration: a method for the prediction of survival outcomes using complex longitudinal and high-dimensional data

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/pencal)](https://cran.r-project.org/package=pencal)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/pencal?color=orange)](http://cranlogs.r-pkg.org/badges/grand-total/pencal)
[![Monthly Downloads](http://cranlogs.r-pkg.org/badges/pencal)](http://cranlogs.r-pkg.org/badges/pencal)

<img src="https://user-images.githubusercontent.com/20061736/162180793-072613f0-a93e-4ef6-b0c4-b8d8a45d770a.png" align="right" alt="" width="250" />

## `pencal`: what's that?
`pencal` is the `R` package that implements Penalized Regression Calibration (PRC), a statistical method that we proposed in Signorelli *et al.* (2021). Penalized regression calibration: A method for the prediction of survival outcomes using complex longitudinal and high-dimensional data. *Statistics in Medicine*, 40 (27), 6178-6196. 
You can read and download the paper (with open access) here: https://onlinelibrary.wiley.com/doi/10.1002/sim.9178

## About this repository
This repository contains the data and code to reproduce the simulations presented in Signorelli et al. (2021). It is divided into 3 subfolders:

* `simulations_PRC_LMM` contains the scripts used for simulations 1 to 6;
* `simulations_PRC_MLPMM` contains the scripts used for simulations 7 to 12;
* `supplementary` contains the results of the simulations presented in the supplementary material.

## How to download pencal
`pencal` is an `R` package that can be downloaded from `CRAN`. This can be done in `R` through the command 
```
install.packages('pencal')
```

If you encounter problems with packages on which `pencal` depends, you may alternatively install the package using `BiocManager`:

```
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install('pencal')
```

## Further information
More information on PRC can be found at the following pages:
* [CRAN package page](https://cran.r-project.org/web/packages/pencal/index.html);
* [my personal website](https://mirkosignorelli.github.io/r.html);
* [vignette that illustrates with examples how to use `pencal`](https://cran.r-project.org/web/packages/pencal/vignettes/vignette.pdf).

A read-only mirror of the package's source code is available at https://github.com/cran/pencal.

