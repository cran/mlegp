`calcLogLike` <-
function(gp) {
	return (dmvnorm(t(gp$Z), matrix(gp$mu, gp$numObs), solve(gp$invVarMatrix), log=TRUE))
}

