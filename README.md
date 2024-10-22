What is HAMT?
======

HAMT is a new Heteroskedasticity Adjusted Multiple Testing procedure that avoids data reduction by standardization, and directly incorporates the side information from
the variances into the testing procedure. Our approach relies on an improved nonparametric empirical Bayes deconvolution estimator that offers a practical strategy for capturing the dependence between the inferential parameter of interest and the variance of the test statistic.

The main idea
------------
Suppose $X_i$, $i=1, \cdots, m$, are independent summary statistics arising from the following random mixture model:

$$X_{i}=\mu_i+\sigma_i\epsilon_{i},~~\epsilon_{i}\stackrel{i.i.d.}{\sim}\eta(\cdot),\quad\mu_i\mid\sigma_i\stackrel{ind.}{\sim}g_\mu(\cdot\mid\sigma_i),\quad\sigma_i\stackrel{i.i.d.}{\sim}g_\sigma(\cdot),$$

where $\epsilon_{i}$ are i.i.d. $0$ mean random variables with a known probability density function (PDF) $\eta(\cdot)$ and $g_\mu(\cdot\mid\sigma_i)$, $g_\sigma(\cdot)$ are, respectively, the PDFs of the unknown mixing distributions of $\mu$ given $\sigma_i$ and $\sigma_i$. Upon observing the pair $(X_i, \sigma_i)$, the goal is to simultaneously test the following $m$ hypotheses: $H_{0, i}: \mu_i \in \mathcal{A_0} \quad \text{versus} \quad  H_{1, i}: \mu_i \notin \mathcal{A_0},~i=1,\ldots,m$, 
where $\mathcal A_0$ represents the indifference region such that the researcher is indifferent to the effects in $\mathcal A_0$. To solve this multiple testing problem involving composite nulls, the proposed HAMT procedure constructs local false discovery rate (Lfdr) statistics by learning $g_\mu(\cdot\mid\sigma_i)$ from the data. The main strategy involves discretizing the support of $g_\mu(\cdot\mid\sigma_i)$ into $S$ equispaced points and then approximating it as a linear combination of $K$ basis functions. The weights $\mathcal W$ in this approximation are determined by solving the following convex optimization problem

$$\min_{\mathcal W}||\mathbf{\varphi}-\mathbf{A}\mathbf{W}||_2^2 \text{ subject to }CW=1,\quad BW\ge 0.$$

Here the matrices \mathbf{\varphi}, A, B and C are known and discussed in Section 3 of the paper. The above problem can be difficult to solve as numerical instabilities may prevent the solver (MOSEK in this case) from satisfying the equality constraints $CW=1$. We use two alternative strategies in our implementation. 

How to use this repository?
----------

This repository holds the scripts that reproduce the analysis in the paper [1]. In particular, the subfolder `simulations` has the R code to reproduce figures 4 to 32 in the paper while the R code for reproducing figures 1 to 3 are available in the main directory. In each case, the code `funcs.R` (available inside the `simulation` subfolder) must be available in the current R working directory. You will also need `MOSEK 9.3 or higher` and the associated R interface available in the [`Rmosek`](https://docs.mosek.com/latest/rmosek/index.html) package.

If the interest is in testing `HAMT` on your own examples then `example 1.R` and `example 2.R` are good places to begin. Again, please make `funcs.R` is in the current working directory. 

References
=======
[1.] [Large-Scale Multiple Testing of Composite Null Hypotheses Under Heteroskedasticity](https://arxiv.org/abs/2306.07362)   
Gang, B. and Banerjee, T.(2024) _(under review)_
