\name{nice}
\alias{nice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ FANOVA formatting functions }
\description{
  Prints out FANOVA decomposition tables in a 'nice' format by rounding entries to the specified number of decimal places
}
\usage{
nice(x, digits = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an FANOVA decomposition table, obtained from \code{FANOVADecomposition} }
  \item{digits}{ integer indicating the precision to be used }
}
\details{
  Rounds elements in all columns of \code{x}, except for the first column, to the specified number of decimal places indicated by \code{precision}
}
\value{
  A \sQuote{nice} version of the object \code{x}  
}
\references{ \url{http://users.nsula.edu/dancikg/mlegp/} }
\author{ Garrett M. Dancik \email{ dancikg@nsula.edu } }
\seealso{ \code{\link{FANOVADecomposition}} for an example }
\keyword{ print }

