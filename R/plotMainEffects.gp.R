`plotMainEffects.gp` <-
function(gp, ylab = "predicted output", graphStyle = 2, verbose = FALSE, no.plot = FALSE, effects = gp$params, length.out = 21, lower = NULL, upper = NULL, FANOVA = FALSE, ...) {

	main = "Main Effects"
	main = ""
	xlab = "param value" # will be overwritten if graphStyle = 3

	includeR2 = TRUE

	if (is.null(lower)) lower = apply(gp$X,2,min)
	if (is.null(upper)) upper = apply(gp$X,2,max)
	effects = toParamIndexes(effects, gp$params)
	param.names = gp$params[effects]

	if (FANOVA) {
		R2 = 100
		if (length(gp$params) > 1) {
			cat("calculating FANOVA...\n")
			FVar = FunctionalVariance(gp, lower = lower, upper = upper)
			R2 = matrix(0, length(effects))
			for (i in 1:length(effects)) {
			   R2[i] = ANOVAMainEffect(gp, effects[i], lower = lower, upper = upper)/FVar * 100
			}
			R2 = format(round(R2, 3))
		}
	}

	index = seq(min(gp$X[,effects]), max(gp$X[,effects]), length.out = length.out)
	preds = matrix(0,ncol = length.out, nrow = length(effects))
	scaleIndex = FALSE
	for (i in 1:length(effects)) {
		if (verbose) {
	        	s = paste("calculating effect for", param.names[i])
	        	cat(s)
			cat("...\n")
		}

		if (min(index) != min(gp$X[,effects[i]]) || max(index) != max(gp$X[,effects[i]])) {
			scaleIndex = TRUE
		} 
       	 	preds[i,] = calcMainEffect(gp, effects[i], length.out = length.out, lower = lower, upper = upper)
  	}

	if (no.plot) {
		if (scaleIndex) index = seq(0,1,length.out=length.out)
		d = list(index = index, preds = preds) 
		return (d)
	}

	if (graphStyle > 2) {
		createWindow(length(param.names))
	}

	if (graphStyle == 2) {
		par(mfrow = c(1,2))
	}

	range = c(min(preds), max(preds))
	for (i in 1:length(effects)) {
	   if (graphStyle == 3) {
		xlab = param.names[i]
	   }
	   index = seq(lower[effects[i]], upper[effects[i]], length.out=length.out)
	   if(scaleIndex) index = seq(0,1,length.out=length.out)
	   xlim = c(min(index),max(index))
	   plot(index, preds[i,], ylab = ylab, type = "n", ylim = range, xlim = xlim, main = main, xlab = xlab)
	   main = ""; xlab = ""; ylab = ""
	   par(xaxt = "n"); par(yaxt = "n") 
           if (graphStyle < 3) par(new=TRUE)
	   lines(index, preds[i,], col=i, lty=i)
           if (graphStyle < 3) par(new=TRUE)
	}

	if (graphStyle == 1 || graphStyle == 2) {
		if (graphStyle == 2) {   ## create the legend plot
		  par(new=FALSE)
		  plot(seq(0,5),seq(0,5), type="n", xlab = "", ylab = "", axes=FALSE)
		}
		cex = 1 #.75
		numEffects = length(effects)
		blank = rep(" ", numEffects)
		temp <- legend("topright", legend = blank,
        		text.width = strwidth("1,000,000"),
		        lty = 1:numEffects, col=1:numEffects, xjust = 1, yjust = 1,
        		title = "Legend",cex = cex )
		if (FANOVA) param.names = paste(paste(param.names, R2, sep=" ("),")", sep="")
		text(temp$rect$left + temp$rect$w, temp$text$y,param.names, pos=2, cex = cex)
		par(new=FALSE)
	}
	par(xaxt = "s"); par(yaxt = "s") 

}

