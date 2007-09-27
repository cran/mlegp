`ANOVAMainEffect` <-
function(gp, effectNum, lower = apply(gp$X, 2, min), upper = apply(gp$X,2,max)) {

  ## calculate cor matrix by integrating out all elements except for effectNum
  corNewObsMatrix = matrix(0, ncol=dim(gp$X)[1])
  for (p in 1:dim(gp$X)[1]) {
    corNewObsMatrix[p] = computeCorElementExcept(gp$beta, gp$a, as.double(gp$X[p,]), effectNum, lower, upper)
  }

  meanRegNoXj = expectedMeanRegExcept(gp, effectNum, lower = lower, upper = upper) 

  integratedValue = integrate(IntegrateMainEffectSquared,lower[effectNum],upper[effectNum],
			gp = gp, meanRegNoXj = meanRegNoXj,
			effectNum = effectNum, corMatrixExceptXj = corNewObsMatrix,  
		        W = gp$invVarMatrix%*%(gp$Z-gp$mu), mu0 = OverallMean(gp, lower = lower, upper = upper))

  return (integratedValue$value * 1 / (upper[effectNum]-lower[effectNum]))

}

