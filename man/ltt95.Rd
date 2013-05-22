\name{ltt95}
\alias{ltt95}
\title{Creates a (1-\eqn{\alpha})-percent CI for a set of LTTs}
\usage{
ltt95(trees, alpha=0.05, log=FALSE, method=c("lineages","times"),
	mode=c("median","mean"), ...)
}
\arguments{
	\item{trees}{is an object of class \code{"multiPhylo"} containing a list of phylogenetic trees.}
	\item{alpha}{confidence level.}
	\item{log}{logical value indicating whether or not to plot on the log-scale.}
	\item{method}{plot the CI on the number of lineages given time (\code{"lineages"}); or on times given a number of lineages (\code{"times"}).}
	\item{mode}{plot the median or mean LTT.}
	\item{...}{optional arguments. So far, \code{res} gives the number of time-steps (defaults to \code{res=100}.}
}
\description{
	This function computes LTT plots for a set of trees & plots a (1-\eqn{\alpha})-percent CI by various methods.
}
\details{
	This function creates a plot and returns a matrix with the results.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ltt}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{diversification analysis}
