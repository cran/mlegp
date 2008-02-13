\name{plotMainEffects}
\alias{plotMainEffects}
\alias{plotMainEffects.gp}
\alias{plotMainEffects.gp.list}
\alias{plotMainEffectsUsingPCWeights}
\alias{calcMainEffect}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sensitivity Analysis for Gaussian Processes and Gaussian Process Lists }
\description{
  Plots main effects of Gaussian process design parameters
}
\usage{
plotMainEffects(gp, ylab = "predicted output", graphStyle = 2, verbose = FALSE, no.plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gp}{ an object of class 'gp' or 'gp.list' }
  \item{ylab}{ optional y-axis label }
  \item{graphStyle}{ integer indicating the format of the graph, as folllows:\tabular{ll}{\tab 0 - no legend\cr
			\tab 1 - legend is on top right of graph\cr
			\tab 2 - legend is to the right of the graph (default)}				
 }
  \item{verbose}{ set \code{verbose} to \sQuote{TRUE} to display status updates while running }
  \item{no.plot}{ set \code{no.plot} to \sQuote{TRUE} to turn plotting off and to return the main effects }
  \item{\dots}{ see \code{Additional argument} sections }
}

\details{

Implements main effects plots described in Schonlau and Welch (2006), using a prior distribution of all components that is ~ independent U(a,b), where (a,b) corresponds to the (min, max) value of that component in the design. For all parameters, (a,b) can be overwritten via the arguments \code{lower} and \code{upper}.

For a single Gaussian process, main effects of the specified parameters are plotted; for Gaussian process lists, the main effects of a single parameter for (a subset of) all Gaussian processes in the list are plotted; for functional output, the main effect is caculated by calculating the main effect for each Gaussian process and premultiplying by the UD matrix to obtain the main effect on the functional output.

Main effects for the parameter \code{x} of non-functional output will plot the main effect of \code{x} evaluated at values seq(min(x), max(x), length.out), where \code{length.out} can be specified by the user. For functional output, \code{plotMainEffects} considers values of \code{x} with \code{length.out} equal to 3; the parameter values can also be specified directly via the argument \code{param.values}.	

If main effects for several parameters are plotted on the same graph, and the ranges for parameters are not identical, they are automatically sclaed to be between 0 and 1.

Additional arguments (described below) will vary depending on whether \code{gp} is of type \code{gp} or type \code{gp.list}.
}

\value{
   If \code{no.plot} is \code{TRUE}, a list with components:
	\item{index}{vector of parameter values}
	\item{preds}{vector of corresponding predicted values}
}

\section{Additional Arguments (General)}{
\item{\code{FANOVA}}{ if TRUE, and graphStyle is 1 or 2, FANOVA results (\% of total functional variance of each effect) is reported in the legend; this option is not implmented for functional output }
\item{\code{lower}}{ a vector for the minimum values of all parameters of the \code{gp} design matrix, defaulting to the minimum value of each parameter in the design }
\item{\code{upper}}{ a vector for the maximum values of all parameters of the \code{gp} design matrix, defaulting to the maximum value of each parameter in the design }
}

\section{Additional Arguments (for single Gaussian processes)}{
\item{\code{effects}}{ a vector of parameter names (or numbers) to plot main effects for, corresponding to the column of the design matrix (default is all parameters) }
\item{\code{length.out}}{ the number of values to plot for each parameter }
}

\section{Additional Arguments (for Gaussian process lists)}{
\item{\code{effects}}{ the name (or number) of the parameter to plot main effects for; only one main effect can be plotted at a time  }
\item{\code{gp.nums}}{ a vector of the Gaussian process numbers to plot (all Gaussian processes in list are analyzed by default) }
\item{\code{length.out}}{ the number of values of each parameter to plot (but does not apply for main effects of functional output }
\item{\code{PC.weights}}{ if the component \code{UD} of \code{gp} is not \code{NULL}, set PC.weights to \code{FALSE} (the default) to plot the main effects of observed response; otherwise the main effects of the principle component weights are plotted } 
\item{\code{param.values}}{ if calculating the main effect of output using principle component weights, the parameter values to look at }
\item{\code{xlab}}{ if calculating the main effect of output using principle component weights, optional x-axis label (default is 'time')  }
\item{\code{xValues}}{ for plotting main effects of functional output, the values to plot along the x-axis }
}




\references{ Schonlau, M. and Welch, W. 2006. Screening the Input Variables to a Computer Model Via Analysis of Variance and Visualization, in Screening: Methods for Experimentation in Industry, Drug Discovery, and Genetics. A Dean and S. Lewis, eds. (New York: Springer). 


\url{http://www.public.iastate.edu/~gdancik/mlegp/}

}
\author{ Garrett M. Dancik \email{ gdancik@iastate.edu } }
\seealso{ \code{\link{FANOVADecomposition}} for quantifying the importance of main effects }
\examples{

## fit a Gaussian process ##
x1 = rep(-5:5,3); x2 = kronecker(1:3, rep(1,11))
y = sin(x1) + 2*x2 + rnorm(length(x1),sd=.001)
fit = mlegp(cbind(x1,x2), y)

## plot all main effects
plotMainEffects(fit) 

\dontrun{   ## library adapt must be loaded for FANOVA=TRUE 
## plot all main effects with FANOVA contributions in legend
plotMainEffects(fit, FANOVA = TRUE) 
}

## plotMainEffects for x2 'manually'
main = plotMainEffects(fit, effects = 2, no.plot=TRUE)
plot(main$index, main$preds, xlab = "x2", ylab = "predicted output")

## fit multiple Gaussian processes ##
x = seq(-5,5)
z1 = 10 - 5*x + rnorm(length(x))
z2 = 4 - 5*x + rnorm(length(x))
z3 = 7*sin(x) + rnorm(length(x))
fitMulti = mlegp(x, cbind(z1,z2,z3), param.names = "x")

# for each gaussian process, plot main effects for parameter 1
plotMainEffects(fitMulti, effects = 1)

# plot the main effect of the first parameter for gaussian processes 1 and 3
plotMainEffects(fitMulti, effects = 1, gp.nums = c(1,3))


## fit multiple Gaussian processes using principle component weights ##
x = seq(-5,5,by=.2)
p = 1:50
y = matrix(0,length(p), length(x))
for (i in p) {
        y[i,] = sin(x) + i + rnorm(length(x), sd  = .01)
}

numPCs = 2
fitPC = mlegp(p, t(y), PC.num = numPCs)

# plot main effects on the functional output ##
plotMainEffects(fitPC, effects = 1, xValues = x)

}

\keyword{ hplot }
\keyword{ htest }

