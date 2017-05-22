## Check consistency of various parts of the implementation.
## Requires that user has implemented SIMULATE code.
##
## 1. what="joint": Report p-value for correct simulation (validates
## simulation code).
## 2. what="marginal": Report percentwise bias of estimates and std
## errors due to Laplace inaccuracy.
##
## Tests do not depend on data. Can be run for any parameter in the
## parameter space. By default use 'last.par.best' if exists -
## otherwise 'par'.

checkConsistency <- function(obj,
                             par = NULL,
                             what = c("marginal", "joint"),
                             hessian = FALSE,
                             n = 10,
                             reDoCholesky = TRUE
                             ) {
    ## Args to construct copy of 'obj'
    args <- as.list(obj$env)[intersect(names(formals(MakeADFun)), ls(obj$env))]
    ## Determine parameter and full parameter to use
    r0 <- r <- obj$env$random
    if( is.null(par) ) {
        ## Default case: Optimization has been carried out by user
        if (is.null(obj$env$last.par.best)) {
            stop("'par' not specified.")
        }
        parfull <- obj$env$last.par.best
        if( any(r) ) par <- parfull[-r] else par <- parfull
    } else {
        ## Custom case: User specifies parameter vector (fixed effects)
        parfull <- obj$env$par
        if( any(r) ) parfull[-r] <- par else parfull <- par
    }
    ## Get names of random effects (excluding profiled parameters)
    if(any(obj$env$profile)) {
        r0 <- r[ ! as.logical(obj$env$profile) ]
        names.profile <- unique(names(parfull[r[as.logical(obj$env$profile)]]))
    } else {
        names.profile <- NULL
    }
    names.random <- unique(names(parfull[r0]))
    ## Use 'parfull' for new object
    args$parameters <- obj$env$parList(par, par = parfull)
    ## Fix all profiled parameters
    map.profile <- lapply(args$parameters[names.profile], function(x)factor(x*NA))
    args$map <- c(args$map, map.profile)
    ## Find randomeffects character
    args$random <- names.random
    args$regexp <- FALSE
    ## Create new object
    newobj <- do.call("MakeADFun", args)
    doSim <- function(...) {
        newobj$env$data <- newobj$simulate(newobj$env$par, complete=TRUE)
        ## Check that random effects have been simulated
        haveRandomSim <- all( names.random %in% names(newobj$env$data) )
        if (haveRandomSim) {
            newobj$env$parameters[names.random] <- newobj$env$data[names.random]
        }
        if(reDoCholesky)
            newobj$env$L.created.by.newton <- NULL
        newobj$env$retape()
        ans <- list()
        if (haveRandomSim) {
            ans$gradientJoint <- newobj$env$f(order=1)
        }
        ans$gradient <- newobj$gr(par)
        if (hessian) ans$hessian <- optimHess(par, newobj$fn, newobj$gr)
        ans
    }
    ans <- lapply(seq_len(n), doSim)
    attr(ans, "par") <- par
    class(ans) <- "checkConsistency"
    ans
}

print.checkConsistency <- function(x, alpha=.05, ...) {
    cat("Parameters used for simulation:\n")
    print(attr(x, "par"))
    ## Check simulation
    check <- function(name = "gradientJoint", get=c("p.value", "bias")) {
        get <- match.arg(get)
        if(is.null(x[[1]][[name]])) {
            return(NA)
        }
        nsim <- length(x)
        mat <- do.call("cbind", lapply(x, function(x)as.vector(x[[name]])))
        mu <- rowMeans(mat)
        n <- length(mu)
        bias <- p.value <- NULL
        if(nsim < n) {
            ## Independence
            warning("Assuming independence in simulation test ", nsim,"<",n)
            q <- apply(mat, 1, function(x)mean(x)/sd(x))
            p.value <- 1 - pchisq(q, df=1)
            p.value <- mean(p.value)
            bias <- apply(mat, 1, function(x)mean(x)/(x))
        } else {
            ## Variance of score = Information
            H <- var(t(mat))
            iH <- try(solve(H), silent=TRUE)
            if(is(iH, "try-error")) {
                warning("Failed to invert information matrix")
                bias <- attr(x, "par") * NA
                p.value <- NA
            } else {
                q <- as.vector( t(mu) %*% iH %*% mu )
                p.value <- 1 - pchisq(q, df=nrow(H))
                bias <- -iH %*% mu
            }
        }
        bias <- as.vector(bias)
        names(bias) <- names(attr(x, "par"))
        list(p.value=p.value, bias=bias)[[get]]
    }
    p.value <- check("gradientJoint", "p.value")
    cat("\n")
    cat("Test correct simulation (p.value):\n")
    print(p.value)
    sim.ok <- p.value > alpha
    if(is.na(sim.ok))
        cat("Full simulation was not available\n")
    else if(!sim.ok)
        cat("Simulation does *not* appear to be correct !!!\n")
    else
        cat("Simulation appears to be correct\n")
    ## Check Laplace:
    cat("\n")
    cat("Estimated parameter bias due to Laplace approximation:\n")
    bias <- check("gradient", "bias")
    print(bias)
}

if(FALSE) {
    library(TMB)
    runExample("sam", exfolder="../../tmb_examples")
    set.seed(123)
    qw <- checkConsistency(obj, opt$par, n=100)
    print.checkConsistency(qw)
    runExample("ar1_4D", exfolder="../../tmb_examples")
    set.seed(123)
    qw <- checkConsistency(obj, opt$par, n=100)
    qw
}
