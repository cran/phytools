\name{treeSlice}
\alias{treeSlice}
\title{Slices the tree at a particular point and returns all subtrees}
\usage{
treeSlice(tree, slice, trivial=FALSE)
}
\arguments{
  \item{tree}{is a phylogenetic tree in \code{"phylo"} format.}
  \item{slice}{a real number indicating the height above the root at which to slice the tree.}
  \item{trivial}{a logical value indicating whether or not to return subtrees with a number of tips less than two (default is \code{FALSE}).}
}
\description{
  This function slices a tree at a particular height above the root and returns all subtrees or all non-trivial subtrees (i.e., subtrees with more than 1 taxon). Uses \code{\link{extract.clade}} in the "ape" package.
}
\value{
  a \code{"phylo"} or \code{"multiPhylo"} object.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{extract.clade}}
}
\keyword{phylogenetics}
\keyword{utilities}
