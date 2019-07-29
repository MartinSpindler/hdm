# hdm: High-Dimensional Metrics
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/hdm)](https://cran.r-project.org/package=hdm)

Implementation of selected high-dimensional statistical and econometric methods for estimation and inference. Efficient estimators and uniformly valid confidence intervals for various low-dimensional causal/ structural parameters are provided which appear in high-dimensional approximately sparse models. Including functions for fitting heteroscedastic robust Lasso regressions with non-Gaussian errors and for instrumental variable (IV) and treatment effect estimation in a high-dimensional setting. Moreover, the methods enable valid post-selection inference and rely on a theoretically grounded, data-driven choice of the penalty. 


## Getting Started with `hdm`

R is an open source software project and can be freely downloaded from the CRAN website along with its associated documentation. There are two options to install the `R` package `hdm` - either installation of the development version or the stable release available at [CRAN](https://CRAN.R-project.org/package=hdm). 

### Development Version

The current development version of the `hdm` package is maintained in this repository and can be installed by the command `devtools::install_github("MartinSpindler/hdm")`. Note that the `devtools` package is required for this command.

### Stable Release

The stable package release is available at [CRAN](https://CRAN.R-project.org/package=hdm). The stable release version can be installed by typing `install.packages("hdm")` in `R`.

### Getting Started: Vignette

After installation, users can get started by following [the package vignette](doc/hdm.pdf).
 
 
## References

V. Chernozhukov, C. Hansen and M. Spindler (2016). "hdm: High-dimensional metrics." arXiv preprint arXiv:1608.00354 (2016), [available online](https://arxiv.org/abs/1603.01700).

A. Belloni, D. Chen, V. Chernozhukov and C. Hansen (2012). Sparse models and methods for optimal instruments with an application to eminent domain. Econometrica 80 (6), 2369-2429.

A. Belloni, V. Chernozhukov and C. Hansen (2013). Inference for high-dimensional sparse econometric models. In Advances in Economics and Econometrics: 10th World Congress, Vol. 3: Econometrics, Cambridge University Press: Cambridge, 245-295.

A. Belloni, V. Chernozhukov, C. Hansen (2014). Inference on treatment effects after selection among high-dimensional controls. The Review of Economic Studies 81(2), 608-650.
