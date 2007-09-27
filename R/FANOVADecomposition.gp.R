`FANOVADecomposition.gp` <-
function(gp, Interaction=TRUE, verbose = TRUE, outputs = NULL, 
		lower = apply(gp$X,2,min), upper = apply(gp$X,2,max), maxpts = NULL) {

	if (verbose) cat("calculating total functional variance...\n")
	FVar = FunctionalVariance(gp, lower = lower, upper = upper)

	numParams = gp$numDim
	numEffects = numParams
	params = gp$params

	if (Interaction) {
		numInts = choose(numParams,2)
		numEffects = numParams + numInts
	}

	R2 = matrix(0,numEffects)
	names = matrix(0,numEffects)
	names[1:numParams] = params

        for (i in 1:numParams) {
	        if (verbose) {
			s = paste("calculating effect #", i)
			s = paste(s, " / ")
			s = paste(s, numEffects)
			cat(s)
			cat("\n")
		}
                R2[i] = ANOVAMainEffect(gp, i, lower = lower, upper = upper) / FVar * 100.0
        }

	if (!Interaction) return (data.frame(param=names, "% contribution"=R2, check.names=FALSE))

        XInts = matrix(0,choose(numParams,2))
        count = numParams + 1
        int1 = 1
        int2 = 2
        for (int1 in 1:(numParams-1)) {
 	       for (int2 in (int1+1):numParams) {
		  names[count] = paste(params[int1], params[int2], sep=":")
			if (verbose) {
				s = paste("calculating effect #", count)
				s = paste(s, " / ")
				s = paste(s, numEffects)
                		cat(s)
				cat("\n")
			}
                 	R2[count] = ANOVAInteractionEffect(gp, c(int1,int2), maxpts = maxpts) / FVar * 100.0
	                count = count + 1
                }
          }

	return (data.frame(param=names, "% contribution"=R2, check.names=FALSE))
}

