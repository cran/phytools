\name{splitTree}
\alias{splitTree}
\title{Split tree at a point}
\usage{
splitTree(tree, split)
}
\arguments{
  \item{tree}{phylogenetic tree.}
  \item{split}{split encoded as a list with two elements: \code{node}: the node number tipward of the split; and \code{bp}: the position along the branch to break the tree, measured from the rootward end of the edge.}
}
\description{
  Primarily an internal function for \code{\link{posterior.evolrate}}, this function splits the tree at a given point, and returns the two subtrees as an object of class \code{"multiPhylo"}.
}
\details{
  Probably do not use this unless you can figure out what you are doing.
}
\value{
  two trees in a list.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
