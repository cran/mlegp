`FunctionalVariance` <-
function(gp, lower = apply(gp$X,2,min), upper = apply(gp$X,2,max) ,mu0 = NULL) {
  if (is.null(mu0)) {
          mu0 = OverallMean(gp, lower = lower, upper = upper)
  }
  numParams = dim(gp$X)[2]
  W = gp$invVarMatrix%*%(gp$Z-gp$mu)
  h = adapt(numParams, lower = lower, upper = upper,functn = IntegrateForVariance, gp = gp, 
	W = W, mu0 = mu0) 
  return (h$value * 1 / prod(upper-lower))
}

