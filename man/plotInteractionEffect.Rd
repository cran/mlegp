\name{plotInteractionEffect}
\alias{plotInteractionEffect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sensitivity Analysis for Gaussian Processes and Gaussian Process Lists }
\description{
  Displays a contour plot that shows the 2-way interaction effect of two parameters from a Gaussian process
}
\usage{
plotInteractionEffect(gp, effects, length1.out = 21, length2.out = 21, 
	lower = apply(gp$X, 2, min), upper = apply(gp$X, 2, max), no.plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gp}{ an object of class \code{gp} }
  \item{effects}{ a vector containing 2 elements corresponding to the parameter numbers or names to plot }
  \item{length1.out}{ the number of values to consider for the first parameter }
  \item{length2.out}{ the number of values to consider for the second parameter }
  \item{lower}{ a vector of minimum values for ALL parameters of the \code{gp} design matrix }
  \item{upper}{ a vector of maximum values for ALL parameters of the \code{gp} design matrix }
  \item{no.plot}{ set \code{no.plot} to \code{TRUE} to turn plotting off and to return the interaction effects }
}
\details{

   An interaction effect of two parameters, x and y, is the predicted output when x and y are fixed, and the remaining parameters are integrated out according to a prior distribution that is independent U(a,b) for all parameters, where for each parameter, (a,b) corresponds to the (min,max) of that parameter in the design. These values can also be overwritten by the arguments \code{lower} and \code{upper}.  

   This function produces a contour plot of predictions for all combinations of the parameters in \code{effects}, where the values of the first parameter are 
   seq(lower[effects[1]], upper[effects[1]], length.out = length1.out) and values for the second parameter are
   seq(lower[effects[2]], upper[effects[2]], length.out = length2.out). 

   For elements not in \code{effects}, \code{lower} and \code{upper} will determine the limits of integration.   

}
\value{
   If \code{no.plot} is \code{TRUE}, a list with components:
	\item{index1}{vector of parameter values for the first parameter }
	\item{index2}{vector of parameter values for the second parameter}
	\item{preds}{a matrix containing the predicted values to be plotted by \code{contour}}
}
\references{ Schonlau, M. and Welch, W. 2006. Screening the Input Variables to a Computer Model Via Analysis of Variance and Visualization, in Screening: Methods for Experimentation in Industry, Drug Discovery, and Genetics. A Dean and S. Lewis, eds. (New York: Springer).  

\url{http://www.public.iastate.edu/~gdancik/mlegp/}
}
\author{ Garrett M. Dancik \email{ dancikg@nsula.edu } }
\seealso{ \code{\link{FANOVADecomposition}} for quantifying the importance of an interaction effect
}
\examples{
## fit the Gaussian process ##
x1 = kronecker(seq(0,1,by=.25), rep(1,5))
x2 = rep(seq(0,1,by=.25),5)
z = 4 * x1 - 2*x2 + x1 * x2 + rnorm(length(x1), sd = 0.001)
fit = mlegp(cbind(x1,x2), z, param.names = c("x1", "x2"))

## plot the interaction effect ##
plotInteractionEffect(fit, effects = c(1,2))

## plot the interaction effect 'manually' ##
int = plotInteractionEffect(fit, effects = c(1,2), no.plot = TRUE)
contour(int$index1, int$index2, int$preds, xlab = "x1", ylab = "x2")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ htest }
