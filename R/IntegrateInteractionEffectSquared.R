`IntegrateInteractionEffectSquared` <-
function(x, gp, meanRegNoX1, meanRegNoX2, meanRegNoX1X2, 
		corMatrixNoX1, corMatrixNoX2, corMatrixNoX1X2, W, mu0, effectNums) {
	if (gp$constantMean == 1) {
		y1 = meanRegNoX1 + gp$sig2 * 
		  (corMatrixNoX1 * 
	 	        exp(-gp$beta[effectNums[1]]*abs(x[1] - gp$X[,effectNums[1]])**gp$a[effectNums[1]])) %*%W
		y2 = meanRegNoX2 + gp$sig2 * 
	  	  (corMatrixNoX2 * 
			exp(-gp$beta[effectNums[2]]*abs(x[2] - gp$X[,effectNums[2]])**gp$a[effectNums[2]])) %*%W
		y12 = meanRegNoX1X2 + gp$sig2 *
		(corMatrixNoX1X2 * 
			exp(-gp$beta[effectNums[1]]*abs(x[1] - gp$X[,effectNums[1]])**gp$a[effectNums[1]]) *
			exp(-gp$beta[effectNums[2]]*abs(x[2] - gp$X[,effectNums[2]])**gp$a[effectNums[2]]) ) %*% W
	}
	else {
		y1 = meanRegNoX1 + x[1]*gp$Bhat[effectNums[1]+1] + gp$sig2 * 
	  	  (corMatrixNoX1 * 
			exp(-gp$beta[effectNums[1]]*abs(x[1] - gp$X[,effectNums[1]])**gp$a[effectNums[1]])) %*%W
		y2 = meanRegNoX2 + x[2]*gp$Bhat[effectNums[2]+1] + gp$sig2 * 
	  	  (corMatrixNoX2 * 
			exp(-gp$beta[effectNums[2]]*abs(x[2] - gp$X[,effectNums[2]])**gp$a[effectNums[2]])) %*%W
		y12 = meanRegNoX1X2 + x[1]*gp$Bhat[effectNums[1]+1] + x[2]*gp$Bhat[effectNums[2]+1] + gp$sig2 *
		  (corMatrixNoX1X2 * 
			exp(-gp$beta[effectNums[1]]*abs(x[1] - gp$X[,effectNums[1]])**gp$a[effectNums[1]]) *
			exp(-gp$beta[effectNums[2]]*abs(x[2] - gp$X[,effectNums[2]])**gp$a[effectNums[2]]) ) %*% W
	}
	return ((y12 - y1 - y2 + mu0)**2)
}

