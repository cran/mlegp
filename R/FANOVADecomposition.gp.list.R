`FANOVADecomposition.gp.list` <-
function(gp, Interaction=TRUE, verbose = TRUE, outputs=1:gp$numGPs, maxpts = NULL, lower = NULL, upper = NULL) {
	if (min(outputs) < 1 || max(outputs) > gp$numGPs) {
		stop("output values must be between 1 and ", gp$numGPs)
	}

	ans = NULL
	for (i in outputs) {
		if (verbose) cat("GP #", i, "\n")
		if (is.null(lower)) lower = apply(gp[[i]]$X, 2, min)
		if (is.null(upper)) upper = apply(gp[[i]]$X, 2, max)
		nextAns = FANOVADecomposition(gp[[i]], Interaction = Interaction, verbose=verbose, lower = lower, upper = upper)
		if (is.null(ans)) {
			ans = nextAns
		}
		else {
			ans = cbind(ans, matrix(nextAns[,2]))
		} 	
	}

	colnames(ans) = c("params",gp$names[outputs])
	return (ans)
}

