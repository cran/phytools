\name{drop.tip.simmap}
\alias{drop.tip.simmap}
\title{Drop tip from SIMMAP format tree}
\usage{
drop.tip.simmap(tree, tip)
}
\arguments{
  \item{tree}{a phylogenetic tree in modified \code{"phylo"} format (see \code{\link{read.simmap}}).}
  \item{tip}{name or names of species to be dropped.}
}
\description{
  This function drops one or multiple tips from the modified \code{"phylo"} object with a mapped binary or multistate trait (see \code{\link{read.simmap}}) while maintaining the matrix \code{$mapped.edge} and list of mappings by branch \code{maps}.
}
\value{
  a modified \code{"phylo"} object containing the elements \code{maps} and \code{$mapped.edge} with the time spent in each state along each edge of the tree.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{brownie.lite}}, \code{\link{drop.tip}}, \code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{sim.history}}
}
\keyword{phylogenetics}
\keyword{utilities}
