\name{threshState}
\alias{threshState}
\title{Computes value for a threshold character from a liability and thresholds}
\usage{
	threshState(x, thresholds)
}
\arguments{
	\item{x}{liability.}
	\item{thresholds}{a named vector containing the thresholds.}
}
\description{
	Primarily to be used internally by \code{\link{ancThresh}}; can also be used to simulate threshold traits.
}
\value{
	A discrete character value.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ancThresh}}, \code{\link{threshDIC}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{simulation}
