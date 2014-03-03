\name{rateshift}
\alias{rateshift}
\title{Find the temporal position of one or more rate shifts}
\usage{
rateshift(tree, x, nrates=1, niter=10, ...)
}
\arguments{
 	\item{tree}{object of class \code{"phylo"}.}
	\item{x}{vector of phenotypic trait values for species. \code{names(x)} should contain the species names and match \code{tree$tip.label}.}
	\item{nrates}{number of rates.}
	\item{niter}{number of iterations of optimization routine to ensure convergence.}
	\item{...}{optional arguments.}
}
\description{
	Function finds the location of one or more rate shifts.
}
\value{
	An object of class \code{"rateshift"}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownie.lite}}
}
\keyword{phylogenetics}
\keyword{comparative method}
