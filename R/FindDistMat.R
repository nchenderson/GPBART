FindDistMat <- function(X, weight.vec=NULL) {
  ## X is the n x p "design" matrix
  ## Can we prevent cases where d_{ij} = 1?
  nvars <- ncol(X)
  DistMat <- matrix(0.0, nrow=nrow(X), ncol=nrow(X))
  if(is.null(weight.vec)) {
    weight.vec <- rep(1/ncol(X), ncol(X))
  }
  for(u in 1:nvars) {
    Amax <- outer(X[,u], X[,u], FUN="pmax")
    Amin <- outer(X[,u], X[,u], FUN="pmin")
    rr <- FindCutPoints(X[,u])
    Gtmp <- stepfun(rr$qq, c(0,cumsum(rr$probs)), right=TRUE)
    
    Gdiff <- matrix(Gtmp(Amax), nrow=nrow(X)) - matrix(Gtmp(Amin), nrow=nrow(X))
    print(Gdiff[1,2])
    DistMat <- DistMat + weight.vec[u]*Gdiff
  }
  return(DistMat)
}
