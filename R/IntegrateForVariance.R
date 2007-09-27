`IntegrateForVariance` <-
function(x,gp, W, mu0) {
        corMatrix = calcCorOneObs(gp$X, gp$beta, gp$a, x)
	if (gp$constantMean == 1) {
		return (   (gp$mu[1] + gp$sig2*corMatrix%*%W - mu0)**2)        
	}
	return (  (c(1,x)%*%gp$Bhat + gp$sig2*corMatrix%*%W - mu0)**2)
}

