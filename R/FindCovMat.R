FindCovMat <- function(X, hyper.pars, weights=NULL) {
  # Compute DistMatrix  with FindDistMat function first
  
  alpha <- hyper.pars[1]
  bbeta <- hyper.pars[2]
  sigmasq.mu <- hyper.pars[3]
  DistMat <- FindDistMat(X, weights=weights) 

  onemD <- 1 - DistMat ## maybe need to consider cases where dij = 1
  AA <- sigmasq.mu/onemD
  BB <- matrix(1.0, nrow=nrow(X), ncol=nrow(X))
  for(k in 1:10) {
    tmp <- ((alpha^k)/(factorial(k)^bbeta))*(onemD^k)
    BB <- BB + tmp
  }
  prior.cov.mat <- AA*(1 - DistMat*BB)
  return(prior.cov.mat)
}