`plotInteractionEffect` <-
function(gp, effects, length1.out = 21, length2.out = 21, 
	lower = apply(gp$X,2,min), upper = apply(gp$X,2,max), no.plot = FALSE) {

	if (length(effects) != 2) {
		stop("effects must be of length 2 for interaction effect")
	}

	if (gp$numDim <= 1) {
		stop("cannot calculate an interaction effect when GP has only 1 parameter!")
	}

  	if(!suppressWarnings(require(adapt, quietly=TRUE, warn.conflicts=FALSE))) {
		stop("package 'adapt' required for plotInteractionEffect")
	}

	index1 = seq(lower[effects[1]], upper[effects[1]], length.out=length1.out)
	index2 = seq(lower[effects[2]], upper[effects[2]], length.out=length2.out)

	N1 = length(index1)
	N2 = length(index2)

	preds = matrix(0,nrow=N1, ncol=N2)
	newTheta = rep(0, dim(gp$X)[2])
	for (i in 1:N1) {
		for (j in 1:N2) {
			corNewObsMatrix = matrix(0, ncol=dim(gp$X)[1])
			newTheta[effects[1]] = index1[i]
			newTheta[effects[2]] = index2[j]
			## compute correlation matrix by integrating out other terms
			for (p in 1:dim(gp$X)[1]) {
				corNewObsMatrix[p] = computeCorElement(gp$beta, gp$a, as.double(gp$X[p,]), newTheta, 
					effects, lower, upper)
			}
	        	preds[i,j] = expectedMeanReg(gp, effects, newTheta[effects], lower = lower, upper = upper) + 
				gp$sig2 * corNewObsMatrix %*% gp$invVar %*% (gp$Z - gp$mu)
		}
	}
	
	if (no.plot) {
		return (list(index1 = index1, index2 = index2, preds = preds))
	}

	contour(index1, index2, preds, xlab = gp$params[effects[1]], ylab = gp$params[effects[2]])
}

