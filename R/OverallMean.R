`OverallMean` <-
function(gp, lower = apply(gp$X,2,min), upper = apply(gp$X,2,max)) {
        ## integrate over ALL terms ##
	numObs = dim(gp$X)[1]
        holdIndex = dim(gp$X)[2]+2
        corNewObsMatrix = matrix(0, ncol=numObs)
        for (p in 1:numObs) {
                corNewObsMatrix[p] = computeCorElement(gp$beta, gp$a, as.double(gp$X[p,]), NULL, holdIndex, lower, upper)
        }
        pred = expectedMeanReg(gp, lower = lower, upper=upper) + 
		gp$sig2 * corNewObsMatrix %*% gp$invVarMatrix %*% (gp$Z - gp$mu)
        return (pred)
}

