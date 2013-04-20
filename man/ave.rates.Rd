\name{ave.rates}
\alias{ave.rates}
\title{Average the posterior rates}
\usage{
	ave.rates(tree, shift, tips, sig1, sig2, ave.shift, showTree=TRUE)
}
\arguments{
	\item{tree}{a tree.}
	\item{shift}{the shift point for this sample.}
	\item{tips}{tip names tipward of \code{shift}.}
	\item{sig1}{rate 1.}
	\item{sig2}{rate 2.}
	\item{ave.shift}{average shift from all samples.}
	\item{showTree}{logical value indicating whether to plot the rate-stretched tree.}

}
\description{
	Primarily internal function for \code{\link{posterior.evolrate}}.
}
\value{
	a list of the rates.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{evol.rate.mcmc}}, \code{\link{minSplit}}, \code{\link{posterior.evolrate}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{bayesian}