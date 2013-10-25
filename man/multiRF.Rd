\name{multiRF}
\alias{multiRF}
\title{Computes Robinson-Foulds distance between a set of trees}
\usage{
multiRF(trees)
}
\arguments{
	\item{trees}{object of class \code{"multiPhylo"} consisting of two or more fully bifurcating, unrooted trees. If trees are rooted, they will be unrooted.}
}
\description{
	Computes the Robinson-Foulds distance between a set of trees in an object of class \code{"multiPhylo"}. If \code{trees} contains a large number of phylogenies (say 100 or 1,000) this will be much faster than calling \code{RF.dist} in the phangorn package for all pairwise comparisons because it avoids repeating some internal calculations. Nonetheless for large numbers \code{multiRF} is slow, and will use lots of memory.
}
\details{
	Computes the Robinson-Foulds distance between all phylogenies in an object of class \code{"multiPhylo"}. Uses \code{\link{prop.part}} internally for most of the heavy lifting.
}
\value{
	A distance matrix.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
