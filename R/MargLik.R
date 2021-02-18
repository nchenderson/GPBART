MargLik <- function(pars, y, X, nu=3) {
    alpha <- pars[1]
    bbeta <- pars[2]
    sigma.mu <- pars[3]
    lambda <- pars[4]
    
    # Compute Covariance matrix bSigma here
    # This will use values of (alpha, beta, sigma.mu)
    #SS <- bSigma
    integrand.fn <- function(u) {
        diag(SS) <- diag(bSigma) + u   
        ldet <- determinant(SS, logarithm=TRUE)$modulus
        ssterm <- crossprod(y, solve(SS, y))
        lsigprior <- (nu + 2)*log(u) + (nu*lambda)/u
        ret <- exp(-ldet/2 - ssterm/2 - lsig.prior/2)
        return(ret)
    }
    ans <- integrate(integrand.fn, lower=0, upper=Inf)$value
    return(ans)
}