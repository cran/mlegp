`plotMainEffectsUsingPCWeights` <-
function(gp, ylab = "predicted output", graphStyle = 2, verbose = FALSE, no.plot = FALSE, holdIndex, param.values = NULL, 
xlab = "time", lower = NULL, upper = NULL, xValues = 1:dim(gp$UD)[1], ...) {
	
	main = "Main effect of "
	main = paste(main, gp$params[holdIndex])
	#main = ""

	if (is.null(param.values)) param.values = seq(min(gp[[1]]$X[,holdIndex]), max(gp[[1]]$X[,holdIndex]),length.out=3)

	UD = gp$UD

        N = length(param.values)
        w = matrix(0,nrow = N, ncol = gp$numGPs)     # GASP output corresponding to weights of SVD
        preds = matrix(0,N)
        output = matrix(0, ncol = dim(UD)[1],nrow = N)

        for (j in 1:gp$numGPs) { ## we need to get weights for all principle components ###
		if (verbose) {
			s = paste("calculating main effect for PC weight #", j)
			cat(s)
			cat("\n")
		}
		if (is.null(lower)) lower = apply(gp[[i]]$X, 2, min)
		if (is.null(upper)) upper = apply(gp[[i]]$X, 2, max)
		w[,j] = calcMainEffect(gp[[j]], holdIndex, index = param.values, lower = lower, upper = upper)
		
        } ## end weights for PC

	for (j in 1:N) {
        	output[j,] = gp$UD%*%w[j,]
	}

	if (no.plot) return (list(index = param.values, preds = output))

	if (graphStyle == 2) {
		par(mfrow = c(1,2))
	}

        for (i in 1:dim(output)[1]) {
                plot(xValues, output[i,], ylim = c(min(output), max(output)), xlab = xlab, ylab = ylab, type = "n", main=main)
		main = ""; xlab = ""; ylab = ""
		par(xaxt = "n"); par(yaxt = "n")
                par(new=TRUE)
        }

        for (i in 1:dim(output)[1]) {
                lines(xValues, output[i,], lty = i)
                par(new=TRUE)
        }
        par(new=FALSE)
	if (graphStyle == 2) {
		par(new=FALSE)
		plot(seq(0,5),seq(0,5), type="n", xlab = "", ylab = "", axes=FALSE)
	}
	
	if (graphStyle > 0) {	
		blank = rep(" ", N)
		temp <- legend("topright", legend = blank,
        		text.width = strwidth("1,000,000"),
		        lty = 1:N, xjust = 1, yjust = 1,
	        	title = "Legend (p1)")
		text(temp$rect$left + temp$rect$w, temp$text$y,param.values, pos=2)
	}
	par(xaxt = "s"); par(yaxt = "s")
}

