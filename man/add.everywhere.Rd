\name{add.everywhere}
\alias{add.everywhere}
\title{Add tip to all edges in a tree}
\usage{
add.everywhere(tree,tip.name)
}
\arguments{
  \item{tree}{an object of class \code{"phylo"}.}
  \item{tip.name}{a string containing the name of the tip to add.}
}
\description{
  This function adds a tip to all branches to the tree and returns a
  a list of trees as an object of class \code{"multiPhylo"}.
}
\details{
}
\value{
  a list of trees as an object of class \code{"multiPhylo"}.  Since the 
  tip can be added to any branch, the length of the list is equal to the
  number of edges in \code{tree}.
}
\references{
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{allFurcTrees}}, \code{\link{exhaustiveMP}}
}
\examples{
}
\keyword{phylogenetics}
\keyword{utilities}
