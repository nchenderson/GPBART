FindCutPoints <- function(x, tau.max=100) {
  ## Function for finding the potential splitting values
  ## for a given covariate column vector x
  
  w <- sort(unique(x))
  nw <- length(w)
  if(nw > tau.max) {
    pvec <- seq(1/tau.max, (tau.max - 1)/tau.max, by=1/tau.max)
    qq <- c(min(x), quantile(x, prob=pvec))
    probs <- 1/tau.max
    ## need to change this a bit
  } else if( nw <= tau.max) {
    qq <- w[-nw]
    probs <- 1/(nw - 1)
  }
  return(list(qq=qq, probs=probs))
}