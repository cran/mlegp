`expectedMeanReg` <-
function (gp, holdIndex = gp$numDim+2, val = NULL, 
	lower = apply(gp$X,2,min), upper = apply(gp$X,2,max)) {

	if (gp$constantMean == 1) {
		return (gp$mu[1])
	}
	## meanReg = XB
	coeffs <- gp$Bhat[2:(gp$numDim+1)]
	coeffs = coeffs / 2 * (upper+lower)
	if (holdIndex[1] > gp$numDim) {              ## overall meanReg
		return (gp$Bhat[1] + sum(coeffs) )
	}
	coeffs[holdIndex] = coeffs[holdIndex] * 2 / (upper+lower)[holdIndex] * val   ## main effect mean reg
	return (gp$Bhat[1] + sum(coeffs))
}

