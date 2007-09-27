`expectedMeanRegExcept` <-
function(gp, exceptIndex, lower = apply(gp$X,2,min), upper = apply(gp$X,2,max)) {
	return (expectedMeanReg(gp, holdIndex = exceptIndex, val = rep(0,length(exceptIndex)), lower, upper )  )
}

