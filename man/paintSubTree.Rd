\name{paintSubTree}
\alias{paintSubTree}
\title{Paint sub-trees with a discrete character}
\usage{
paintSubTree(tree,node,state,anc.state="1",stem=FALSE)
}
\arguments{
  \item{tree}{a phylogenetic tree as an object of class \code{"phylo"} or a modified object with mapped character traits.}
  \item{node}{an integer specifying the node number tipward of which the function should paint the derived state.}
  \item{state}{a string (or numeric value) specifying the state to paint on the tree tipward of \code{node}.}
  \item{anc.state}{the ancestral state to use; will only be applied if there are presently no character values mapped on the tree.}
  \item{stem}{logical or numeric value indicating whether to use the derived state on the stem leading to \code{node} (or not, if \code{stem=FALSE}), or, alternatively, what fraction of the stem should be assigned to the derived clade. Note that for tip clades \code{stem=FALSE} is not allowed.}
}
\description{
	This function maps an arbitrary discrete character history on the tree.
}
\value{
  a modified phylogenetic tree of class \code{"phylo"} with the following additional elements:
  \item{maps}{a list of named vectors containing the times spent in each state on each branch, in the order in which they occur.}
  \item{mapped.edge}{a matrix containing the total time spent in each state along each edge of the tree.}
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{plotSimmap}}, \code{\link{sim.history}}
}
\keyword{phylogenetics}
\keyword{comparative method}
