\name{minRotate}
\alias{minRotate}
\title{Rotates all nodes of the tree to minimize the difference in order with a vector}
\usage{
minRotate(tree, x)
}
\arguments{
	\item{tree}{tree.}
	\item{x}{numeric vector.}
}
\description{
	This function rotates all the nodes of the tree to try and minimize the different between the order of the tips and the rank-order of a numeric vector \code{x}.
}
\details{
	Primarily designed to be used internally by \code{\link{phylo.to.map}}.
}
\value{
	A object of class \code{"phylo"}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
