`mlegp` <-
function(X, Z, constantMean = 1, nugget = NULL, min.nugget = 0, param.names = NULL, gp.names = NULL, 
	PC.UD = NULL, PC.num = NULL, PC.percent = NULL, 
	simplex.ntries = 5, simplex.maxiter = 100, simplex.abstol = 1e-16, simplex.reltol = 1e-8,  
	BFGS.maxiter = 500, BFGS.tol = 0.01, BFGS.h = 1e-10, seed = 0, verbose = 1) {

	X = as.matrix(X)
	Z = as.matrix(Z)

	if (!is.null(param.names) && length(param.names) != dim(X)[2]) {
		stop("length of param.names must match number of columns of X")
	} 

	if (!is.null(PC.percent)) PC.num = numSingularValues(Z, PC.percent)

	if (!is.null(PC.num)) {
		weights = pcweights(Z,weights.num = PC.num)
		Z = t(weights$Vprime)
		if (PC.num == 1) {
			stop("error: currently not implemented for PC.num = 1 (set PC.num directly)")
		} 
		PC.UD = weights$UD
	}

	if (!is.null(dim(Z))){
		numGPs = dim(Z)[2]
	}
	else {
		numGPs = 1
		Z = matrix(Z)
	}

	if (is.null(gp.names)) {
		gp.names = paste("gp #", 1:numGPs)
	}

	if (length(gp.names) != numGPs) {
		stop("length of gp.names must match number of GPs")
	}

	if (dim(X)[1] != dim(Z)[1]) {
		stop("error: X matrix and output matrix must have same number of rows; \n\tfor PC weights, make sure each column of Z is a computer model run")
	}

	if (min(diag(var(Z))) <= 0) {
		stop("error: Z does not vary for at least one column")
	}

	if (dim(Z)[2] > 1 && !is.null(nugget) && length(nugget) > 1) {
		stop("error: cannot handle nugget matrix with multiple field observations; try fitting 1 gp at a time")
	}

	if (!is.null(nugget) && length(nugget) > dim(Z)[1]) {
		stop("length of nugget matrix must match number of observations")
	}

	if ( min(apply(X, 2, max) - apply(X, 2, min)) == 0) {
		stop("error: at least one parameter does not vary in design")
	}

	if (constantMean != 1) {
		ones = rep(1,dim(X)[1])
		dX = cbind(ones, X)
		t = try(solve(t(dX)%*%dX),TRUE)
		if (class(t) == "try-error") {
			stop("error: design matrix with intercept is not full rank; set constantMean = 1")
		}
	}


	if (!is.null(nugget) && length(nugget) == 1 && nugget <= 0) {
		## check for ANY reps
		if (anyReps(X)) stop("at least 2 or more inputs are identical...must use nugget!")
	}
	
	if (is.null(nugget)) {
		if (!anyReps(X)) {
			#nugget = 0
			if (verbose > 0) cat("no reps detected - nugget will not be estimated\n")
		}	
		else {
			if (verbose > 0) cat("reps detected - nugget will be estimated\n")	
		}
	}	

	numEstimates = dim(X)[2] + 1  ## correlation parameters + sig2GP 
	if (anyReps(X) || !is.null(nugget) || length(nugget) > 1) 
		numEstimates = numEstimates + 1   ## add nugget term if needed
	if (constantMean == 1) {
		numEstimates = numEstimates + 1
	}
	else {
		numEstimates = numEstimates + dim(X)[2] + 1
	}

	estimates = rep(0,numEstimates)

	if (verbose > 0) cat("\n\n")
	
	l = NULL
	for (i in 1:numGPs) {
	    if (verbose > 0) {
		    cat("========== FITTING GP #")
		    cat(i)
		    cat("==============================\n")
	    }

	    if (is.null(nugget) || (length(nugget) == 1 && nugget != 0)) {
		if (anyReps(X)) nugget = estimateNugget(X,Z[,i] )	
		
	    }

	    success = 0
	    estimates = rep(0,numEstimates)
 	    returnFromC = .C("fitGPfromR", as.double(X), as.integer(nrow(X)), as.integer(ncol(X)),
		as.double(Z[,i]), as.integer(nrow(Z)), 
		as.integer(constantMean), 
		as.integer(simplex.ntries), as.integer(simplex.maxiter), 
			as.double(simplex.abstol), as.double(simplex.reltol),
		as.integer(BFGS.maxiter), as.double(BFGS.tol), as.double(BFGS.h), 
		as.integer(seed), as.double(nugget), as.integer(length(nugget)), as.double(min.nugget),
		estimates = as.double(estimates), verbose = as.integer(verbose), success = as.integer(success), 
		PACKAGE="mlegp")

	    if (returnFromC$success != 0) {
		cat("ERROR: GP cannot be created\n")
	        return (NULL)
            }	
	    estimates = returnFromC$estimates

    	    numParams = dim(X)[2]
	    regSize = 1
	    meanReg = estimates[1]
	    if (constantMean == 0) {
		regSize = dim(X)[2] + 1
		meanReg = estimates[1:regSize]
	    }

 	    ## remove meanReg params; now we have correlation params, sig2, and (possibly) nugget
	    estimates = estimates[(regSize+1):length(estimates)]   ## remove mean reg params

	    beta = estimates[1:numParams]
	    a = rep(2, numParams)
	    sig2 = estimates[numParams+1]

	    if (is.null(nugget)) nugget = 0

	    if (length(nugget) > 1) {
		if (verbose > 0) {
			cat(paste("nugget matrix will be scaled by: ", estimates[numParams+2]))
			cat("\n")
		}
		nugget = nugget * estimates[numParams+2]
	    }
	    else {
		if (nugget > 0) {
			nugget = estimates[numParams+2]
		}
	    }
	
	   nugget = nugget + min.nugget

 	    if (is.null(param.names)) param.names = paste("p",1:dim(X)[2],sep="")

	    if (verbose > 0) cat("creating gp object...")
	    fit1 = createGP(X, as.matrix(Z[,i]), beta, a, meanReg, sig2, nugget, 
		param.names = param.names, constantMean = constantMean)
	    if (verbose > 0) cat("...done\n\n")

	    l[[i]] = fit1

	    }

	if (numGPs == 1) {
		return (l[[1]])
	}
	return (gp.list(l, UD = PC.UD, param.names = param.names, gp.names = gp.names))
}

