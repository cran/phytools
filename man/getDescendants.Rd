\name{getDescendants}
\alias{getDescendants}
\title{Get descendant node numbers}
\usage{
getDescendants(tree,node,curr=NULL)
}
\arguments{
  \item{tree}{a phylogenetic tree as an object of class \code{"phylo"}.}
  \item{node}{an integer specifying a node number in the tree.}
  \item{curr}{the set of previously stored node numbers - used in recursive function calls.}
}
\description{
	This function returns the set of node & tip numbers descended from \code{node}.
}
\value{
  a vector.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{paintSubTree}}
}
\keyword{phylogenetics}
\keyword{utilities}
