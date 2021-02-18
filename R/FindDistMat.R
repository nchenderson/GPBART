FindDistMat <- function(X, tau.max=100, weights=NULL) {
  ## X is the n x p "design" matrix
  ## Can we prevent cases where d_{ij} = 1?
  nvars <- ncol(X)
  DistMat <- matrix(0.0, nrow=nrow(X), ncol=nrow(X))
  if(is.null(weights)) {
    weights <- rep(1/ncol(X), ncol(X))
  }
  for(u in 1:nvars) {
    Amax <- outer(X[,u], X[,u], FUN="pmax")
    Amin <- outer(X[,u], X[,u], FUN="pmin")
    rr <- FindCutPoints(X[,u], tau.max=tau.max)## add tau.max
    Gtmp <- stepfun(rr$qq, c(0,cumsum(rr$probs)), right=TRUE)
    
    Gdiff <- matrix(Gtmp(Amax), nrow=nrow(X)) - matrix(Gtmp(Amin), nrow=nrow(X))
    DistMat <- DistMat + weights[u]*Gdiff
  }
  return(DistMat)
}
