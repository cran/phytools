\name{drop.tip.simmap}
\alias{drop.tip.simmap}
\title{Drop tip from SIMMAP format tree}
\usage{
drop.tip.simmap(tree,tip)
}
\arguments{
  \item{tree}{a phylogenetic tree in modified \code{"phylo"} format (see \code{\link{read.simmap}}).}
  \item{tip}{name of species to be dropped.}
}
\description{
  This function drops a tip from the modified \code{"phylo"} object with a mapped binary or multistate trait (see \code{\link{read.simmap}}) while maintaining the matrix \code{$mapped.edge}.
}
\details{
  Note that this function does not, as yet, modify the element \code{maps}; only \code{mapped.edge}.
}
\value{
  a modified \code{"phylo"} object containing the element \code{$mapped.edge} with the time spent in each state along each edge of the tree.
}
\references{
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{brownie.lite}}, \code{\link{drop.tip}}, \code{\link{evol.vcv}}, \code{\link{read.simmap}}
}
\examples{
}
\keyword{phylogenetics}
