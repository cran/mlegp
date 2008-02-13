\name{FANOVADecomposition}
\alias{FANOVADecomposition}
\alias{FANOVADecomposition.gp}
\alias{FANOVADecomposition.gp.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sensitivity Analysis for Gaussian Processes and Gaussian Process Lists }
\description{
  Performs Functional Analysis of Variance (FANOVA) decomposition on a Gaussian process or (a subset of) all Gaussian processes in a Gaussian process list 
}
\usage{
FANOVADecomposition(gp, Interaction = TRUE, verbose = TRUE, outputs = NULL, maxpts = NULL, 
	lower = NULL, upper = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gp}{ an object of class \code{gp} or \code{gp.list} }
  \item{Interaction}{ if TRUE (the default), calculates two-way factor interactions}
  \item{verbose}{ if TRUE (the default), prints status updates }
  \item{outputs}{ if gp is of class \code{gp.list}, a vector of integers specifying the Gaussian processes to analyze; all Gaussian processes in the list are analyzed by default }
  \item{maxpts}{ the maximum number of function evaluations used in the calculation of each predicted output, see \code{adapt} }
  \item{lower}{ a vector of minimum values of ALL parameters of the \code{gp} design matrix, defaulting to the minimum value of each parameter in the design }
  \item{upper}{ a vector of maximum values of ALL parameters of the \code{gp} design matrix, defaulting to the maximum value of each parameter in the design }
}
\details{
  Implements the FANOVA decomposition of Schonlau and Welch (2006), for main and two-way factor effects, using a prior distribution of all components that is ~ independent U(a,b), where (a,b) corresponds to the (min,max) value of that
  component in the design. For all parameters, (a,b) can be overwritten via the arguments \code{lower} and \code{upper}.  

  The functions \code{integrate} and \code{adapt} from the library \code{adapt} are used for one dimensional and multidimensional integration, respectively.
 
}
\value{
  a data.frame where the first column contains the names of the effects, and elements in row i and column j correspond to the percentage of the total functional variance of Gaussian process (j-1) accounted for by effect i 
}
\references{ Schonlau, M. and Welch, W. 2006. Screening the Input Variables to a Computer Model Via Analysis of Variance and Visualization, in Screening: Methods for Experimentation in Industry, Drug Discovery, and Genetics. A Dean and S. Lewis, eds. (New York: Springer).  

\url{http://www.public.iastate.edu/~gdancik/mlegp/}
}


\author{ Garrett M. Dancik \email{ gdancik@iastate.edu } }
\note{ 
     This function requires the library \code{adapt}  

     The total contributions from all main and two-way factor interactions should never exceed 100\%; however, it is possible that this will happen due to numerical rounding error   		
 
     If R starts to hang, changing the default value of \code{maxpts} is needed (though convergence no longer guaranteed)

     From \pkg{adapt}: 

     Setting \code{maxpts} to NULL asks the function to keep doubling maxpts
     (starting at \code{max(minpts,500, r(ndim))}) until the desired
     precision is achieved or R runs out of memory...
}
\seealso{ \code{\link{FANOVAFunctions}}; \code{\link{nice}} for rounding entries in the FANOVA table; \code{\link{plotMainEffects}} and \code{\link{plotInteractionEffect}} for plotting main and 2-way interaction effects;
\code{\link{integrate}} and package \pkg{adapt} for details on single and multi-dimensional integration
	  
}
\examples{
\dontrun{ ## Do not run the following code unless the library adapt is loaded ##

## fit the Gaussian process ##
x1 = kronecker(seq(0,1,by=.25), rep(1,5))
x2 = rep(seq(0,1,by=.25),5)
z = 4 * x1 - 2*x2 + x1 * x2 + rnorm(length(x1), sd = 0.001)
fit = mlegp(cbind(x1,x2), z, param.names = c("x1", "x2"))

## Find percent contribution of all effects ##
FANOVADecomposition(fit)

## Round contributions to 3 digits ##
nice(FANOVADecomposition(fit), 3)
} ### end don't run
}
\keyword{ htest }
