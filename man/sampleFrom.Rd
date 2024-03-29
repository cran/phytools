\name{sampleFrom}
\alias{sampleFrom}
\title{Sample from a set of distributions}
\usage{
sampleFrom(xbar=0, xvar=1, n=1, randn=NULL, type="norm")
}
\arguments{
	\item{xbar}{a named vector of means.}
	\item{xvar}{a named vector of variances.}
	\item{n}{a vector containing the sample sizes of each species.}
	\item{randn}{a range of sample sizes are to be random.}
	\item{type}{\code{"norm"} is the only distribution implemented so far.}
}
\description{
	Samples from a set of normal distributions with parameters given in \code{xbar} and \code{xvar}.
}
\value{
	A vector, with labels.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{statistics}
\keyword{utilities}
