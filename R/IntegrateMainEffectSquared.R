`IntegrateMainEffectSquared` <-
function(x, gp, meanRegNoXj, effectNum, corMatrixExceptXj, W, mu0) {
	y = x
	for (i in 1:length(x)) {
		if (gp$constantMean == 1) {
			y[i] = (meanRegNoXj + gp$sig2 * 
	  		  (corMatrixExceptXj * 
				exp(-gp$beta[effectNum]*abs(x[i] - gp$X[,effectNum])**gp$a[effectNum])) %*%W
			  - mu0)**2
		}
		else {
			y[i] = (meanRegNoXj + x[i]*gp$Bhat[effectNum+1] + gp$sig2 * 
	  		  (corMatrixExceptXj * 
				exp(-gp$beta[effectNum]*abs(x[i] - gp$X[,effectNum])**gp$a[effectNum])) %*%W
			  - mu0)**2
		}
	}
	return (y)
}

