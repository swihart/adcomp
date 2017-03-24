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
    class(ans) <- "checkConsistency"
    ans
}

if(FALSE) {
    library(TMB)
    ##runExample("sam", exfolder="../../tmb_examples")
    runExample("ar1_4D", exfolder="../../tmb_examples")
    qw <- checkConsistency(obj, opt$par, n=10)
    mat <- sapply(qw, function(x)x$gradient)
    rowMeans(mat)
}
