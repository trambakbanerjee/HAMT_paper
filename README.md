What is HAMT?
======

HAMT is a new Heteroskedasticity Adjusted Multiple Testing procedure that avoids data reduction by standardization, and directly incorporates the side information from
the variances into the testing procedure. Our approach relies on an improved nonparametric empirical Bayes deconvolution estimator that offers a practical strategy for capturing the dependence between the inferential parameter of interest and the variance of the test statistic.

The main idea
------------

How to use this repository?
----------

This repository holds the scripts that reproduce the analysis in the paper [1]. In particular, the subfolder `simulations` has the R code to reproduce figures 4 to 32 in the paper while the R code for reproducing figures 1 to 3 are available in the main directory. In each case, the code `funcs.R` (available inside the `simulation` subfolder) must be available in the current R working directory. You will also need `MOSEK 9.3 or higher` and the associated R interface available in the [`Rmosek`](https://docs.mosek.com/latest/rmosek/index.html) package.

If the interest is in testing `HAMT` on your own examples then `example 1.R` and `example 2.R` are good places to begin. Again, please make `funcs.R` is in the current working directory. 

References
=======
[1.] [Large-Scale Multiple Testing of Composite Null Hypotheses Under Heteroskedasticity](https://arxiv.org/abs/2306.07362)   
Gang, B. and Banerjee, T.(2024) _(under review)_
