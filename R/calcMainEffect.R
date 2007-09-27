`calcMainEffect` <-
function(gp, holdIndex, index = NULL, length.out=21, lower = apply(gp$X,2,min), upper = apply(gp$X,2,max)) {

	if (is.null(index)) index = seq(lower[holdIndex], upper[holdIndex], length.out=length.out)
        N = length(index)

        preds = matrix(0,N)
        newTheta = rep(0, gp$numDim)

        for (i in 1:N) {
                corNewObsMatrix = matrix(0, ncol=dim(gp$X)[1])
                newTheta[holdIndex] = index[i]
                ## compute correlation matrix by integrating out other terms
                for (p in 1:dim(gp$X)[1]) {
                        #print(p)
                        corNewObsMatrix[p] = computeCorElement(gp$beta, gp$a, as.double(gp$X[p,]), newTheta, 
				holdIndex, lower, upper)

                }
	        preds[i] = expectedMeanReg(gp, holdIndex, newTheta[holdIndex], lower = lower, upper = upper) + 
				gp$sig2 * corNewObsMatrix %*% gp$invVarMatrix %*% (gp$Z - gp$mu)

        }
        return (t(preds))
}

