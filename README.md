# rpl: Randomized pairwise likelihood method for complex statistical inferences

The *rpl* package can be installed by using the following commands:

```
library(devtools)
devtools::install_github("andreamrau/rpl")
library(rpl)
```
To also build the vignette (still under construction), you can use the following (note that this 
will require the installation of some extra packages):

```
devtools::install_github("andreamrau/rpl", build_vignettes=TRUE)
library(rpl)
```

*rpl* implements an algorithm for a randomized pairwise likelihood approach that enables computationally
efficient statistical inference, including the construction of confidence intervals, in cases where 
the full likelihood is too complex to be used (e.g., multivariate count data). The primary functions
of this package are as follows:

- `rpl_optim`: optimize model parameters for a given correlation matrix structure (`"unstructured"`,
`"exchangeable"`, `"block_exchangeable"`, or `"factor"`)
- `rpl_se`: calculate standard errors for model parameters for a given correlation matrix structure (`"unstructured"`)
- `init`: obtain reasonable values to initialize marginal and correlation parameters for use in the `rpl_optim` function
- `simulate_mvt_poisson`: simulate multivariate Poisson data with a given correlation matrix structure 
(`"unstructured"`,
`"exchangeable"`, `"block_exchangeable"`, or `"factor"`)

If you use *rpl* in your research, please cite our work:

- Mazo, G., Karlis, D., Rau, A. (2021) A randomized pairwise likelihood
method for complex statistical inferences. In revision. [hal-03126621](https://hal.archives-ouvertes.fr/hal-03126621)

