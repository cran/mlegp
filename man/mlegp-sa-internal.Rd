\name{mlegp-sa-internal}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{mlegp-sa-internal}
\alias{Integrate1DCor}
\alias{Integrate1DCorSquared}
\alias{IntegrateForVariance}
\alias{IntegrateMainEffectSquared}
\alias{IntegrateInteractionEffectSquared}
\alias{createIntegratedCorMatrixExcept}
\alias{computeCorElement}
\alias{computeCorElementExcept}
\alias{expectedMeanReg}
\alias{expectedMeanRegExcept}
\title{ Internal Sensitivity Analysis Functions For Gaussian Processes }
\description{
	Internal Sensitivity Analysis Functions For Gaussian Processes
}
\details{
	These functions are called by FANOVADecomposition, plotMainEffects, and plotInteractionEffect,
	and should not be called by the user. Some of these require the library \code{adapt}.
}
\references{ \url{http://users.nsula.edu/dancikg/mlegp/} }
\author{ Garrett M. Dancik \email{ dancikg@nsula.edu } }
\keyword{ internal }
\keyword{ htest }
