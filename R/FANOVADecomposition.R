`FANOVADecomposition` <-
function(gp, Interaction = TRUE, verbose=TRUE, outputs = NULL, maxpts = NULL, lower = NULL, upper = NULL) {
	if (gp$numDim == 1) {
		stop("for 1 dimensional GP, the single parameter accounts for 100% of the total variance")
	}
  	if(!suppressWarnings(require(adapt, quietly=TRUE, warn.conflicts=FALSE))) {
		stop("package 'adapt' required for FANOVADecomposition")
	}

	if (is.gp.list(gp)) {
		if (is.null(outputs)) outputs = 1:gp$numGPs
		if (is.null(lower)) lower = apply(gp[[1]]$X, 2, min)
		if (is.null(upper)) upper = apply(gp[[1]]$X, 2, max)
	}
	else {
		if (is.null(lower)) lower = apply(gp$X, 2, min)
		if (is.null(upper)) upper = apply(gp$X, 2, max)
	}
	NextMethod("FANOVADecomposition", gp, Interaction = Interaction, verbose=verbose, outputs = outputs, lower = lower, upper = upper)	
}

