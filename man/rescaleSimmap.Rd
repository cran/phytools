\name{rescaleSimmap}
\alias{rescaleSimmap}
\title{Rescale SIMMAP style tree}
\usage{
	rescaleSimmap(tree, totalDepth=1.0)
}
\arguments{
	\item{tree}{a phylogenetic tree in modified \code{"phylo"} format with a discrete character mapping (e.g., see \code{\link{read.simmap}} or \code{\link{make.simmap}}).}
	\item{totalDepth}{a new desired value for the total depth (height) of the tree. If the tree is not ultrametric this is the maximum height.}
}
\description{
	This function re-scales a SIMMAP style tree to an arbitrary total height.
}
\details{
	Replaces \code{rescaleTree} in the 'geiger' package for SIMMAP style trees.
}
\value{
  a phylogenetic tree in modified \code{"phylo"} format.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{make.simmap}}, \code{\link{read.simmap}}
}
\keyword{phylogenetics}
\keyword{utilities}
