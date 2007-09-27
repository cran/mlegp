\name{FANOVA}
\alias{FANOVAFunctions}
\alias{OverallMean}
\alias{ANOVAMainEffect}
\alias{ANOVAInteractionEffect}
\alias{FunctionalVariance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sensitivity Analysis for Gaussian Processes and Gaussian Process Lists }
\description{
  Functional Analysis of Variance (FANOVA) Functions for Gaussian Processes
}
\usage{
OverallMean(gp, lower = apply(gp$X, 2, min), upper = apply(gp$X, 2, max))
FunctionalVariance(gp, lower = apply(gp$X, 2, min), upper = apply(gp$X, 2, max), mu0 = NULL)
ANOVAMainEffect(gp, effectNum, lower = apply(gp$X, 2, min), upper = apply(gp$X, 2, max))
ANOVAInteractionEffect(gp, effectNums, mu0 = NULL, lower = apply(gp$X, 2, min), upper = apply(gp$X, 2, max), 
	maxpts = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gp}{ an object of type \code{gp} }
  \item{lower}{ a vector of minimum values for ALL parameters of the \code{gp} design matrix }
  \item{upper}{ a vector of maximum values for ALL parameters of the \code{gp} design matrix }
  \item{effectNum}{ the effect number or name of the desired main effect }
  \item{effectNums}{ a vector of 2 elements corresponding to the effect numbers or names of the desired interaction effect }
  \item{mu0}{ the overall mean, useful for speeding up computations}
  \item{maxpts}{ the maximum number of function evaluations used in the calculation of each predicted output, see \code{adapt} }
}
\details{
  Implements the FANOVA decomposition of Schonlau and Welch (2006), for main and two-way factor effects, using a prior distribution of all components that is ~ independent U(a,b), where (a,b) corresponds to the (min,max) value of that
  component in the design. For all parameters, (a,b) can be overwritten via the arguments \code{lower} and \code{upper}.  

  \code{OverallMean} and \code{ANOVAMainEffect} use the function \code{integrate} for integration; FunctionalVariance and ANOVAInteractionEffect also use the library \code{adapt} for multi-dimensional integration.
  
 
}
\value{
  \code{OverallMean} returns the mean value of the Gaussian process predictor

  \code{FunctionalVariance} returns the overall variance of the Gaussian process predictor

  \code{ANOVAMainEffect} returns the variation due to the single effect specified by \code{effectNum}

  \code{ANOVAInteractionEffect} returns the variation due to the joint two-way effect specified by \code{effectNms}
}
\references{ Schonlau, M. and Welch, W. 2006. Screening the Input Variables to a Computer Model Via Analysis of Variance and Visualization, in Screening: Methods for Experimentation in Industry, Drug Discovery, and Genetics. A Dean and S. Lewis, eds. (New York: Springer).  

\url{http://www.public.iastate.edu/~gdancik/mlegp/}

} 
\author{ Garrett M. Dancik \email{ gdancik@iastate.edu } }
\note{ 

     If R starts to hang, changing the default value of \code{maxpts} is needed (though convergence no longer guaranteed)

     From \pkg{adapt}: 

     Setting \code{maxpts} to NULL asks the function to keep doubling maxpts
     (starting at \code{max(minpts,500, r(ndim))}) until the desired
     precision is achieved or R runs out of memory...

}
\seealso{ The user does not have to directly call these functions, and instead should use \code{\link{FANOVADecomposition}}

\code{\link{integrate}} and package \pkg{adapt} for details on single and multi-dimensional integration
}
\examples{

\dontrun{
## do not run the following code unless the library adapt is loaded ##
## fit the Gaussian process ##
x1 = kronecker(seq(0,1,by=.25), rep(1,5))
x2 = rep(seq(0,1,by=.25),5)
z = 4 * x1 - 2*x2 + x1 * x2 + rnorm(length(x1))
fit = mlegp(cbind(x1,x2), z, param.names = c("x1", "x2"))

## calculate the proportion of total variance accounted for by x1 ##

ANOVAMainEffect(fit, "x1") / FunctionalVariance(fit) * 100
## check this against 'FANOVADecomposition' ##
FANOVADecomposition(fit)
}

}
\keyword{ internal }
\keyword{ htest }

