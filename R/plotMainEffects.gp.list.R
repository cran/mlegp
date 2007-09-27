`plotMainEffects.gp.list` <-
function(gp, ylab = "predicted output", graphStyle = 2, verbose = FALSE, no.plot = FALSE, effects, gp.nums = 1:gp$numGPs, length.out = 21, param.values = NULL, xlab = NULL, PC.weights = FALSE, lower = NULL, upper = NULL, FANOVA = FALSE, ...) {

	if (length(effects) > 1) {
		stop ("must select only 1 effect for a gp")
	}
	
	if (is.null(lower)) lower = apply(gp[[1]]$X, 2, min)
	if (is.null(upper)) upper = apply(gp[[1]]$X, 2, max)

	if (!is.null(gp$UD) && !PC.weights) {

		if (is.null(xlab)) xlab = "time"
		return (plotMainEffectsUsingPCWeights(gp, ylab = ylab, graphStyle = graphStyle, 
			verbose = verbose, no.plot = no.plot, effects, param.values = param.values, xlab = xlab, lower = lower, upper = upper, ...) ) 
	}


	if (min(gp.nums) < 1 || max(gp.nums) > gp$numGPs) {
		s = paste("gp.nums  must be between 1 and gp$numGPs = ", gp$numGPs)
		stop (s)
	} 

	effects = toParamIndexes(effects, gp$params)
	if (FANOVA) {
		R2 = rep(100, length(gp.nums))
		if (length(gp$params) > 1) {
			cat("calculating FANOVA...\n")
			for (i in length(gp.nums)) {
				FVar = FunctionalVariance(gp[[gp.nums[i]]], lower = lower, upper = upper)
				R2[i] = ANOVAMainEffect(gp[[gp.nums[i]]], effects, lower = lower, upper = upper)/FVar * 100
			}
			R2 = format(round(R2, 3))
		}
	}
	if (graphStyle < 3) {
		main = "Main Effect of"
		main = paste(main, gp$params[effects])
		main = paste(main, "for multiple gp's")
	}

	
	index = param.values
	if (is.null(index)) index = seq(lower[effects], upper[effects], length.out = length.out)

	preds = matrix(0, ncol = length.out, nrow = length(gp.nums))
	for (i in 1:length(gp.nums)) {
		if (verbose) {
			cat ("calculating main effect for gp #")
			cat(gp.nums[i])
			cat("...\n")
		}
		preds[i,] = calcMainEffect(gp[[gp.nums[i]]], effects, index = index, lower = lower, upper = upper)
	}

	if (no.plot) {
		d = list(index = index, preds = preds)
		return (d)
	}
	
	m1 = min(preds)
	m2 = max(preds)

	if (graphStyle == 2) par(mfrow = c(1,2))
	if (graphStyle > 2) {
		createWindow(length(gp.nums))
	}

	if (is.null(xlab)) xlab = "param value"

	for (i in 1:length(gp.nums)) {

		if (graphStyle == 3) {		
			main = "Main Effect of "
			main = paste(main, gp$params[effects])
			main = paste(main, " for gp #")
			main = paste(main, gp.nums[i])
		}

		plot(index, preds[i,], type = "l", lty = i, ylim = c(m1, m2), col=i, main=main, 
			xlab = xlab, ylab = ylab)
		if (graphStyle < 3) par(new = TRUE)
	}
	
	if (graphStyle == 1 || graphStyle == 2) {
		if (graphStyle == 2) {
		  par(new=FALSE)
		  plot(seq(0,5),seq(0,5), type="n", xlab = "", ylab = "", axes=FALSE)
		}
	
		blank = rep(" ", length(gp.nums))
		temp <- legend("topright", legend = blank,
        		text.width = strwidth("1,000,000"),
		        lty = 1:length(gp.nums), col=1:length(gp.nums), 
			xjust = 1, yjust = 1,
        		title = "Legend")
		txt = paste("gp: ", gp.nums)
		if (FANOVA) txt = paste(paste(txt, R2, sep=" ("),")", sep="")
		text(temp$rect$left + temp$rect$w, temp$text$y,txt, pos=2)
		par(new=FALSE)
	}

}

