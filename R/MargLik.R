MargLogLik <- function(hyper.pars, y, X, sigest, nu=3, H = 200, weights=NULL) {
    ## The input hyperparameters are (alpha, beta, kappa, q)
    alpha <- hyper.pars[1]
    bbeta <- hyper.pars[2]
    kappa <- hyper.pars[3]
    qq <- hyper.pars[4]
    lambda <- ((sigest*sigest)/nu)*qchisq(1 - qq, df=nu)
    print(hyper.pars)
    # Compute Covariance matrix bSigma here
    # This will use values of (alpha, beta, sigma.mu)
    hp <- c(alpha, bbeta, 1/(4*kappa*kappa*H), lambda)
    n <- length(y)
    tmp <- FindCovMat(X, hyper.pars=hp, weights=weights)
    bSigma <- H*tmp
    SS <- bSigma
    integrand.fn <- function(u) {
        ret <- rep(NA, length(u))
        for(k in seq_len(length(u))) {
           diag(SS) <- diag(bSigma) + u[k]   
           ldet <- determinant(SS, logarithm=TRUE)$modulus + n*log(2*pi)
           ssterm <- as.numeric(crossprod(y, solve(SS, y)))
           lsig.prior <- (nu + 2)*log(u[k]) + (nu*lambda)/u[k]
           ret[k] <- exp(-ldet/2 - ssterm/2 - lsig.prior/2)
        }
        return(ret)
    }
    ans <- integrate(integrand.fn, lower=0, upper=Inf)$value
    return(-log(ans) - (nu/2)*log(lambda)) # return negative log-likelihood
}

#FindHypers <- function(y, X) {
    # First transform outcomes appropriately
#    mu.y <- (max(y) + min(y))/2
#    s.y <- max(y) - min(y)
#    zz <- (y - mu.y)/s.y
    ## compute residual ...
    
#    lb <- c(0.5, 0.5, 0.05, 0)
#    ub <- c(0.99, 4, 20, ?)
#    oo <- optim(init.vals, lower=lb, upper=ub, fn=MargLik, y=zz, X=X)$par
#    return(oo)
#}



