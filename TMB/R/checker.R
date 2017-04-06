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
                             par = obj$par,
                             what = c("marginal", "joint"),
                             hessian = FALSE,
                             n = 10
                             ) {
    r <- obj$env$random
    parfull <- obj$env$par
    randomNames <- unique(names(parfull[r])) ## Fixme: Remove profile parameters
    if( any(r) ) parfull[-r] <- par else parfull[] <- par
    parameters <- obj$env$parList(par = parfull)
    doSim <- function(...) {
        obj$env$parameters <- parameters
        obj$env$data <- obj$simulate(complete=TRUE)
        ## Check that random effects have been simulated
        haveRandomSim <- all( randomNames %in% names(obj$env$data) )
        if (haveRandomSim) {
            obj$env$parameters[randomNames] <- obj$env$data[randomNames]
        }
        obj$env$retape()
        ans <- list()
        if (haveRandomSim) {
            ans$gradientJoint <- obj$env$f(order=1)
        }
        ans$gradient <- obj$gr(par)
        if (hessian) ans$hessian <- optimHess(par, obj$fn, obj$gr)
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
            iH <- solve(H)
            q <- as.vector( t(mu) %*% iH %*% mu )
            p.value <- 1 - pchisq(q, df=nrow(H))
            bias <- iH %*% mu
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
    qw
    runExample("ar1_4D", exfolder="../../tmb_examples")
    set.seed(123)
    qw <- checkConsistency(obj, opt$par, n=100)
    qw
}
