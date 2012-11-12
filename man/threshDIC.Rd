\name{threshDIC}
\alias{threshDIC}
\title{Deviance Information Criterion from the threshold model}
\usage{
	threshDIC(tree, x, mcmc, burnin)
}
\arguments{
	\item{tree}{phylogenetic tree.}
	\item{x}{named vector containing discrete character states.}
	\item{mcmc}{list object returned by \code{\link{ancThresh}}.}
	\item{burnin}{number of generations (not samples) to exclude as burn in.}
}
\description{
	This function computes the Deviance Information Criterion from the MCMC object returned by \code{\link{ancThresh}}.
}
\value{
	a scalar
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ancThresh}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{bayesian}
