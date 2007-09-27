#`plotMainEffects` <-
#function(gp, effects = NULL, gp.nums = NULL, length.out = 21, param.values = NULL,xlab = NULL, ylab = "predicted output", graphStyle = 2, verbose = FALSE, no.plot=FALSE, no.UD = FALSE) {
#	UseMethod("plotMainEffects", gp)
#}

`plotMainEffects` <-
function(gp, ylab = "predicted output", graphStyle = 2, verbose = FALSE, no.plot = FALSE, ...) {
	UseMethod("plotMainEffects", gp)
}
