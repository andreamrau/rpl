#' Randomized pairwise likelihood method for complex statistical inferences
#'
#' Randomized pairwise likelihood approach to enable computationally
#' efficient statistical inference, including the construction of
#' confidence intervals, in cases where the full likelihood is too
#' complex to be used (e.g., multivariate count data).
#'
#' \tabular{ll}{ Package: \tab rpl\cr Type: \tab Package\cr Version:
#' \tab 0.0.1\cr Date: \tab 2021-11-08\cr License: \tab GPL (>=3.3.1)\cr LazyLoad:
#' \tab no\cr }
#'
#' @name rpl-package
#' @aliases rpl-package
#' @docType package
#' @author Gildas Mazo, Dimitris Karlis, Andrea Rau
#'
#' Maintainer: Andrea Rau <\url{andrea.rau@@inrae.fr}>
#' @references
#' Mazo, G., Karlis, D., Rau, A. (2021) A randomized pairwise likelihood
#' method for complex statistical inferences. In revision. hal-03126621.
#'
#' @keywords multivariate
#' @importFrom numDeriv grad
#' @importFrom pbivnorm pbivnorm
#' @importFrom copula rCopula normalCopula P2p
#' @importFrom stats qnorm optimize cor optim ppois qpois rbinom runif
#' @importFrom BiocParallel bpparam bplapply
NULL

#' Wrapper for the 'optim' function (includes Fisher transformation)
#'
#' @param param_init Vector providing initial values for all parameters
#' @param mydata Data matrix
#' @param pi Bernoulli sampling parameter (0 < pi <= 1)
#' @param offsets If needed, a matrix of offsets (optional)
#' @param block_indices List of vectors providing indices for each block
#' for blockwise exchangeable correlation matrices
#' @param seed Seed to use for random number generation
#' @param count_pairs If \code{TRUE}, return the number of observations retained for
#' each pair of variables used in the randomized pairwise likelihood
#' @param parallel If \code{TRUE}, some computations are parallelized with the
#' BiocParallel package
#' @param BPPARAM If \code{parallel=TRUE}, the registry of back-ends for use
#' in parallelized calculations
#' @param method Optimization method to be used (see \code{?optim} for more
#' information)
#'
#' @return Optimized values for all parameters using the randomzed
#' pairwise likelihood.
#'
#' @references
#' Mazo, G., Karlis, D., Rau, A. (2021) A randomized pairwise likelihood
#' method for complex statistical inferences. In revision. hal-03126621.
#' @export
#' @example /inst/examples/rpl-package.R
rpl_optim <- function(param_init, mydata, pi, offsets=NULL,
                      block_indices=NULL,
                      seed=12345, count_pairs=FALSE,
                      parallel=FALSE, BPPARAM=NULL,
                      method = "L-BFGS-B") {
    d <- dim(mydata)[2]
    n <- dim(mydata)[1]
    ## Transform correlation parameters to z using Fisher transformation
    param_init[-c(1:d)] <- fisher_trans(param_init[-c(1:d)])
    o <- optim(param_init,
               rpl_loglike,
               mydata = mydata,
               pi=pi,
               offsets = offsets,
               block_indices = block_indices,
               seed=seed,
               count_pairs=count_pairs,
               parallel=parallel,
               BPPARAM = BPPARAM,
               method = method, #"L-BFGS-B",
               # upper = c(rep(Inf, d), rep(1, d*(d-1)/2)),
               # lower = c(rep(0, d), rep(-1, d*(d-1)/2)),
               control = list(maxit = 30000))
    ## Back transform correlation parameters using inverse Fisher transform
    o$par[-c(1:d)] <- inverse_fisher_trans(o$par[-c(1:d)])
    return(o)
}

#' Fisher transformation applied within function
#'
#' par is a vector of mean parameters followed by Fisher-transformed
#' correlation parameters
#'
#' @param par Vector providing parameter values
#' @param mydata Data matrix
#' @param pi Bernoulli sampling parameter (0 < pi <= 1)
#' @param offsets If needed, a matrix of offsets (optional)
#' @param block_indices List of vectors providing indices for each block
#' for blockwise exchangeable correlation matrices
#' @param seed Seed to use for random number generation
#' @param count_pairs If \code{TRUE}, return the number of observations retained for
#' each pair of variables used in the randomized pairwise likelihood
#' @param parallel If \code{TRUE}, some computations are parallelized with the
#' BiocParallel package
#' @param BPPARAM If \code{parallel=TRUE}, the registry of back-ends for use
#' in parallelized calculations
#'
#' @return Negative log-likelihood for user-provided parameter vector,
#' calculated with the randomized pairwise likelihood approach.
#' If \code{count_pairs = TRUE}, counts of the number of observations used
#' for each pair are also returned.
#'
#' @references
#' Mazo, G., Karlis, D., Rau, A. (2021) A randomized pairwise likelihood
#' method for complex statistical inferences. In revision. hal-03126621.
#' @export
#'
rpl_loglike <- function(par,
                        mydata,
                        pi = 1,
                        offsets = NULL,
                        ## list with indices of blocks
                        block_indices = NULL,
                        seed = NULL,
                        count_pairs = FALSE,
                        parallel = FALSE,
                        BPPARAM=ifelse(parallel,bpparam(), NULL)) {

    if(!is.null(seed)) set.seed(seed)
    prob <- pi
    d <- dim(mydata)[2]
    n <- dim(mydata)[1]
    x <- mydata
    u <- v <- matrix(NA, n, d)
    upp <-  1 - 10^-12

    ## Back-transform to correlations
    par[-c(1:d)] <- inverse_fisher_trans(par[-c(1:d)])

    ## Dispersions
    if(length(par) == (d+1)) {
        ## Exchangeable correlation matrix
        par <- c(par[1:d],rep(par[d+1],d*(d-1)/2))
    } else if(length(par) == (d + d)) {
        ## One-factor correlation matrix: corr(i,j) = theta_i * theta_j
        theta <- par[-(1:d)]
        corrmat <- outer(theta, theta)
        diag(corrmat) <- 1
        par <- c(par[1:d], P2p(corrmat)) # lexicographical order
    } else if(length(par) == (d + (d*d-d)/2)) {
        ## Unstructured correlation matrix
        par <- par
    } else if(!is.null(block_indices)) {
        ## Block exchangeable correlation matrix
        b <- length(block_indices)
        if(length(par[-c(1:d)]) != (b*b+b)/2) {
            stop("# blocks not compatible with # correlations")
        }
        if(!all.equal(unlist(block_indices), 1:d)) {
            stop("Please sort data by blocks")
        }
        theta <- par[-c(1:d)]
        corrmat <- matrix(0, nrow = d, ncol = d)
        index <- 1
        for(b1 in 1:b) {
            for(b2 in b1:b) {
                corrmat[block_indices[[b1]], block_indices[[b2]]] <-
                    theta[index]
                corrmat[block_indices[[b2]], block_indices[[b1]]] <-
                    theta[index]
                index <- index +1
            }
        }
        diag(corrmat) <- 1
        par <- c(par[1:d], P2p(corrmat)) # lexicographical order
    } else {
        stop("No other types of correlation matrices currently supported.")
    }

    theta <- createcor(par[-c(1:d)], d)

    ## Means
    if(is.null(offsets)) {
        m <- par[1:d]
        u <- ppois(x - 1, matrix(rep(m, each=nrow(x)), ncol=ncol(x)))
        v <- ppois(x, matrix(rep(m, each=nrow(x)), ncol=ncol(x)))
    } else {
        ## log(E(counts)) = log(offset) + X\beta => E(counts) = exp(log(offset) + X\beta)
        m <- mydata * 0
        for(i in 1:ncol(mydata)) {
            m[,i] <- exp(par[i] + log(offsets[,i]))
        }
        u <- ppois(x - 1, m)
        v <- ppois(x, m)
    }

    u <- (u == 0) * 10 ^ -12 + ((u <= upp) &
                                    (u > 0)) * u  +
        (u > upp) * 0.9999999999
    v <- (v == 0) * 10 ^ -12 + ((v <= upp) &
                                    (v > 0)) * v  +
        (v > upp) * 0.9999999999
    u[is.na(u)] <- 10 ^ -12
    v[is.na(v)] <- 10 ^ -12

    if(prob == 1) {
        ij <- matrix(0, nrow=d, ncol=d)
        ij[upper.tri(ij)] <- 1
        ij <- which(ij == 1, arr.ind=TRUE)
        if(!parallel) {
            compos <- sum(unlist(lapply(1:((d *d - d)/2),
                                        function(xx, ij, u, v, theta) {
                                            contribution(u[,ij[xx,1]],
                                                         u[,ij[xx,2]],
                                                         v[,ij[xx,1]],
                                                         v[,ij[xx,2]],
                                                         theta[ij[xx,1],
                                                               ij[xx,2]])
                                        }, ij=ij, u=u, v=v, theta=theta)))
        } else {
            compos <- sum(unlist(bplapply(1:((d *d - d)/2),
                                          function(xx, ij, u, v, theta) {
                                              contribution(u[,ij[xx,1]],
                                                           u[,ij[xx,2]],
                                                           v[,ij[xx,1]],
                                                           v[,ij[xx,2]],
                                                           theta[ij[xx,1],
                                                                 ij[xx,2]])
                                          }, ij=ij, u=u, v=v, theta=theta,
                                          BPPARAM=BPPARAM)))
        }
        pairchoice <- matrix(1, nrow = nrow(mydata),
                             ncol=(d*(d-1)/2))
    } else if(pi > 0 & pi < 1) {
        if(parallel) {
            message("Note: parallel computing only currently available for pi=1")
        }
        ## Select pairs of variables for each observation
        pairchoice <- matrix(rbinom((d*(d-1)/2)*nrow(mydata),
                                    prob=prob, size=1),
                             nrow = nrow(mydata), ncol=(d*(d-1)/2))
        tmp_mat <- matrix(0, nrow = d, ncol = d)
        tmp_mat[upper.tri(tmp_mat)] <- 1
        colnames(pairchoice) <- apply(which(tmp_mat == 1, arr.ind=TRUE),
                                      1, paste0, collapse=".")

        compos <- 0
        for(i in 1:ncol(pairchoice)) {
            index1 <-
                as.numeric(unlist(lapply(strsplit(colnames(pairchoice)[i],
                                                  split = ".", fixed = TRUE),
                                         function(x) x[1])))
            index2 <-
                as.numeric(unlist(lapply(strsplit(colnames(pairchoice)[i],
                                                  split = ".", fixed = TRUE),
                                         function(x) x[2])))
            obs_index <- which(pairchoice[,i] == 1)
            ## AR: right behavior if a pair is selected for no variables??
            if(length(obs_index) == 0) next;
            compos <- compos + contribution(u[obs_index, index1],
                                            u[obs_index, index2],
                                            v[obs_index, index1],
                                            v[obs_index, index2],
                                            theta[index1, index2])
        }
    } else stop("pi must take a value between 0 and 1")
    if(!count_pairs) {
        return(compos)
    } else {
        return(list(compos = compos,
                    count_pairs = colSums(pairchoice)))
    }
}

#' Calculate standard errors for parameters estimated with the RPL
#'
#' @param mydata Data matrix
#' @param parameters Vector providing estimated values for all parameters
#' @param pi Bernoulli sampling parameter (0 < pi <= 1)
#' @param offsets If needed, a matrix of offsets (optional)
#' @param verbose If \code{TRUE}, include verbose output
#'
#' @return A named list containing the following:
#' \itemize{
#'    \item se - standard errors for each estimated parameter
#'    \item parameters - user-provided estimated parameter values
#' }
#'
#' @references
#' Mazo, G., Karlis, D., Rau, A. (2021) A randomized pairwise likelihood
#' method for complex statistical inferences. In revision. hal-03126621.
#' @export
#' @example /inst/examples/rpl-package.R
#'
rpl_se <- function(mydata, parameters, pi, offsets=NULL, verbose=FALSE) {

    ##  mydata is a n x d matrix of data
    ##  parameters is a vector of estimated parameters
    ##    with d mean parameters at the begining and the
    ##    correlation parameters after that using
    ##    (1,2),(1,3),...,(1,d),(2,3),...,(2,d),...,(d-1,d)
    ##  offsets is a set of offset, same dimension as mydta

    prob <- pi
    off <- offsets
    d <- ncol(mydata)
    n <- nrow(mydata)
    if(is.null(off)) {
        off <- matrix(1,n,d)
    }
    npar <- d + d*(d-1)/2
    if(length(parameters) != npar)
        stop("SE calculation only currently implemented for unstructured correlation matrices.")
    all2 <- NULL
    sde <- matrix(0, npar, npar)
    xx <- mydata  ### just to check things we have to change it later
    theta <- parameters
    metr <- d
    m2 <- 0

    #### for all pairs
    for(i in 1:(d-1) ) {
        for(j in (i+1):d) {
            m2 <- m2+1
            metr <- metr+1
            if(verbose) print(c(i,j))
            sdetemp <- matrix(0, npar, npar)
            nn <- 0

            ### for all observations, try to estimate with Monte Carlo
            for(k in 1:n) {
                temp <- grad(dmass_offset,
                             c(theta[i], theta[j], theta[metr]),
                             x1=xx[k,i], x2=xx[k,j],
                             off1=off[k,i], off2=off[k,j],
                             method="simple")
                #### in case of problems in the calculation we just
                #### get rid of this
                if((all(!is.na(temp))) & (all(temp!=0))) {
                    sdtemp <- rep(0,npar)
                    sdtemp[c(i,j,metr)] <- temp
                    sdetemp <- sdetemp + sdtemp %*% t(sdtemp)
                    nn <- nn + 1
                }
            }
            sde <- sde + sdetemp/nn
        }
    }
    sde <- solve(sde)/(n*prob)
    se <- diag(sde)^0.5
    return(list(se=se, parameters=theta))
}

#' Simulate multivariate Poisson data with a given correlation
#' matrix structure using Gaussian copulas
#'
#' @param n Number of observations to simulate
#' @param d Dimension of variables to simulate
#' @param param Named list (\code{"margins"}, \code{"disp"})
#' of parameters to be used for simulation. Parameters
#' should first include marginal parameters, then vectorized correlation parameters
#' in lexicographical order. The size of the \code{disp} vector depends
#' on the \code{type} of correlation matrix structure to be simulated.
#' @param type Type of correlation matrix structure to be simulated:
#' \code{unstructured} ((d^2 - d)/2 correlation parameters),
#' \code{exchangeable} (1 correlation parameter),
#' \code{one-factor} (d theta parameters, where rho_{ij} = theta_i * theta_j),
#' \code{block_exchangeable} (for b blocks, (b^2 - b)/2 correlation parameters).
#' @param block_indices For \code{type="block_exchangeable"}, list
#' containing the indices of variables belonging to each of the b blocks.
#'
#' @return
#' \item{data }{Matrix of simulated data}
#' \item{param }{True parameter values used to simulate data}
#' @export
#'
simulate_mvt_poisson <- function(n = 1000, d = 4,
                                 param,
                                 type = "unstructured",
                                 block_indices = NULL) {
    if(type == "exchangeable") {
        m <- param$margins
        u <- rCopula(n, normalCopula(param$disp, dim = d))
        mydata <- matrix(NA, nrow = nrow(u), ncol = ncol(u))
        for(i in 1:d) {
            mydata[,i] <- qpois(u[,i], m[i])
        }
    } else if(type == "unstructured") {
        u <- rCopula(n, normalCopula(param$disp, dim = d, dispstr="un"))
        mydata <- matrix(NA, nrow = nrow(u), ncol = ncol(u))
        for(i in 1:d) {
            mydata[,i] <- qpois(u[,i], param$margins[i])
        }
    } else if(type == "factor") {
        corrmat <- outer(param$disp, param$disp)
        u <- rCopula(n, normalCopula(P2p(corrmat), dim = d, dispstr="un"))
        mydata <- matrix(NA, nrow = nrow(u), ncol = ncol(u))
        for(i in 1:d) {
            mydata[,i] <- qpois(u[,i], param$margins[i])
        }
    } else if(type == "block_exchangeable") {
        if(is.null(block_indices)) {
            stop("Block indices needed for block exchangeable sims")
        }
        ## Block exchangeable correlation matrix
        b <- length(block_indices)
        theta <- param$disp
        corrmat <- matrix(0, nrow = d, ncol = d)
        index <- 1
        for(b1 in 1:b) {
            for(b2 in b1:b) {
                corrmat[block_indices[[b1]], block_indices[[b2]]] <-
                    theta[index]
                corrmat[block_indices[[b2]], block_indices[[b1]]] <-
                    theta[index]
                index <- index +1
            }
        }
        u <- rCopula(n, normalCopula(P2p(corrmat), dim = d, dispstr="un"))
        mydata <- matrix(NA, nrow = nrow(u), ncol = ncol(u))
        for(i in 1:d) {
            mydata[,i] <- qpois(u[,i], param$margins[i])
        }
    } else {
        stop("No other simulation strategies currently supported.")
    }
    return(list(data = mydata, param = param))
}


#' Various parameter initialization methods for use with the rpl_optim function
#'
#' @param mydata Data matrix
#' @param type Type of intialization to perform, according to the desired
#' structure of correlation matrix:
#' \code{unstructured} ((d^2 - d)/2 correlation parameters estimated using
#' the empirical Pearson correlation matrix),
#' \code{exchangeable} (1 correlation parameter estimated by averaging off-diagonal
#' of the empirical Pearson correlation matrix),
#' \code{one-factor} (d theta parameters, where rho_{ij} = theta_i * theta_j,
#' estimated using Nelson-Meader optimization of fit compared to the estimated
#' Pearson correlation matrix),
#' \code{block_exchangeable} (for b blocks, (b^2 - b)/2 correlation parameters
#' estimated by averaging empirical Pearson correlation matrix within each block).
#' In all cases, marginal parameters are initalized using the plug-in method
#' with empirical marginal means.
#' Several alternative initializations are also proposed for unstructured
#' correlation matrices: PPearson correlation (\code{unstructured_A},
#' equivalent to \code{unstructured}), pairwise
#' optimization with all data (\code{unstructured_B}), pairwise
#' optimization with random subsample of data (\code{unstructured_C}),
#' random initialization for the copula parameters.
#' (\code{unstructured_D})
#'
#' @param pi Bernoulli sampling parameter (0 < pi <= 1)
#' @param block_indices For \code{type="block_exchangeable"}, list
#' containing the indices of variables belonging to each of the b blocks.
#'
#' @return Vector of parameter values to be used to initialize the
#' \code{rpl_optim} function.
#' @export
#'
init <- function(mydata, type = "unstructured", pi=NULL,
                 block_indices = NULL) {
    n <- nrow(mydata)
    d <- ncol(mydata)
    ## Estimate  Poisson parameters
    PoissonMeans <- colMeans(mydata)
    if(type %in% c("unstructured", "unstructured_A")) {
        param_init <-
            c(PoissonMeans,
              cor(mydata, method="pearson")[lower.tri(cor(mydata,
                                                          method="pearson"))])
    } else if(type == "exchangeable") {
        param_init <-
            c(PoissonMeans,
              mean(cor(mydata, method="pearson")[lower.tri(cor(mydata,
                                                               method="pearson"))]))
    } else if(type == "unstructured_B") {
        ## Plug-in Poisson params + estimate copula params using ALL data
        estimCopCorMat <- matrix(1, ncol=d, nrow=d)
        for (i in 1:(d - 1)) {
            for (j in (i + 1):d) {
                estimCopCorMat[i,j] <- estimCopCorMat[j,i] <-
                    optimize(funToOptim1, interval = c(0,1),
                             PoissonMeans, mydata, i, j)$minimum
            }
        }
        param_init <- c(PoissonMeans, estimCopCorMat[lower.tri(estimCopCorMat)])
    } else if(type == "unstructured_C") {
        ## plug-in Poisson params + estimate copula params using SOME data
        if(is.null(pi)) stop("pi must be provided for type C init")
        estimCopCorMat <- matrix(1, ncol=d, nrow=d)
        for (i in 1:(d - 1)) {
            for (j in (i + 1):d) {
                estimCopCorMat[i,j] <- estimCopCorMat[j,i] <-
                    optimize(funToOptim2, interval=c(0,1),
                             PoissonMeans, mydata, i, j, n, pi)$minimum
            }
        }

        param_init <-  c(PoissonMeans,
                         estimCopCorMat[lower.tri(estimCopCorMat)])
    } else if(type == "unstructured_D") {
        param_init <- c(PoissonMeans, runif(d*(d-1)/2, 0, 1))
    } else if(type == "block_exchangeable") {
        if(is.null(block_indices)) {
            stop("block_indices needed for block_exchangeable initialization")
        }
        tmp <- cor(mydata, method="pearson")
        b <- length(block_indices)
        if(!all.equal(unlist(block_indices), 1:d)) {
            stop("Please sort data by blocks")
        }
        theta <- c()
        index <- 1
        for(b1 in 1:b) {
            for(b2 in b1:b) {
                tmp2 <- tmp[block_indices[[b1]],
                            block_indices[[b2]]]
                theta[index] <- mean(tmp2[upper.tri(tmp2)])
                index <- index +1
            }
        }
        param_init <- c(PoissonMeans, theta)
    } else if(type == "factor") {
        tmpmat <- cor(mydata, method="pearson")
        initinit <- runif(d, min = 0, max = 1)
        o <- optim(
            # par=fisher_trans(rep(0, length=d)),
            par=fisher_trans(initinit),
            fn=fn_factor, tmpmat=tmpmat, method = "L-BFGS-B",
            control=list(maxit = 50000))
        param_init <- c(PoissonMeans, inverse_fisher_trans(o$par))
    } else {
        stop("No other initialization methods currently supported.")
    }
    return(param_init)
}

## NOT EXPORTED: apply inverse Fisher transformation
inverse_fisher_trans <- function(z) {
    return((exp(2 * z) - 1) / (exp(2 * z) + 1))
}

## NOT EXPORTED: apply Fisher transformation
fisher_trans <- function(r) {
    return(1/2 * log((1+r)/(1-r)))
}

## NOT EXPORTED: estimate Gaussian copula
gaussian_copula <- function(u, v, a) {
    ### to eliminate overflows we added the two lines
    u <- (u==0)*(0.000000001) + (u>0)*u
    v <- (v==0)*(0.000000001) + (v>0)*v
    temp <- pbivnorm(qnorm(u), qnorm(v), a)
    t <- !(is.finite(temp))
    return(temp)
}

## NOT EXPORTED: calculate density mass
dmass <- function(u, v, u2, v2, theta) {
    temp <- gaussian_copula(u2, v2, theta) -
        gaussian_copula(u, v2, theta) -
        gaussian_copula(u2, v, theta) +
        gaussian_copula(u, v, theta)
    temp <- (temp <= 0) * (10^-180) + (temp > 0) * temp
    return(temp)
}

## NOT EXPORTED: calculate density mass with offsets (same as dmass otherwise)
dmass_offset <- function(theta,x1,x2,off1,off2) {
    m <- theta[1:2]
    alpha <- theta[3]
    u <- ppois(x1, m[1]*off1)
    v2 <- ppois(x2-1, m[2]*off2)
    u2 <- ppois(x1-1, m[1]*off1)
    v <- ppois(x2, m[2]*off2)
    temp <- gaussian_copula(u2, v2, alpha) -
        gaussian_copula(u, v2, alpha) -
        gaussian_copula(u2, v, alpha) +
        gaussian_copula(u, v, alpha)
    temp <- (temp <= 0) * (10^-180) + (temp > 0) * temp
    return(log(temp))
}


## NOT EXPORTED: calculate negative log-likelihood
contribution <- function(u, v, u2, v2, theta) {
    temp1 <- dmass(u, v, u2, v2, theta)
    logl <- sum(log(temp1))
    return(-logl)
}

## NOT EXPORTED: create a (symmetric) correlation matrix from a
## parameter vector (lexicographic order)
createcor <- function(vec, d, optimized = FALSE) {
    if(!optimized) {
        dd <- d - 1
        mat <- matrix(NA, dd, dd)
        kato <- 1
        metr <- dd
        for (i in 1:dd) {
            mat[i, ] <- c(rep(0, i - 1), vec[kato:metr])
            kato <- metr + 1
            metr <- metr + (dd - i)
        }
        cor <- rbind(cbind(0, mat), 0)
    } else {
        cor <- matrix(0, d, d)
        cor[lower.tri(cor)] <- vec
        cor <- t(cor)
    }
    return(cor)
}

## NOT EXPORTED: internal function for init
funToOptim1 <- function(rho, PoissonMeans, mydata, i, j) {
    rpl_loglike(par=c(PoissonMeans[c(i,j)], rho),
                mydata=mydata[,c(i,j)], pi=1)
}

## NOT EXPORTED: internal function for init
funToOptim2 <- function(rho, PoissonMeans, mydata, i, j, n, pi){
    rpl_loglike(par=c(PoissonMeans[c(i,j)], rho),
                mydata=mydata[sample(1:n,n*pi),c(i,j)],
                pi=1)
}

## NOT EXPORTED: internal function for init
fn_factor <- function(par, tmpmat) {
    par2 <- inverse_fisher_trans(par)
    cor2 <- outer(par2, par2)
    return(sum((tmpmat[lower.tri(tmpmat)] -
                    cor2[lower.tri(cor2)])^2))
}