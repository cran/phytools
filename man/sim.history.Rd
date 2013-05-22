\name{sim.history}
\alias{sim.history}
\title{Simulate stochastic character history under some model}
\usage{
sim.history(tree, Q, anc=NULL, nsim=1)
}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"}.}
	\item{Q}{a matrix containing the instantaneous transition rates between states.}
	\item{anc}{an optional value for the state at the root node; if \code{NULL}, a random state will be assigned.}
	\item{nsim}{number of simulations.}
}
\description{
	This function simulates a stochastic character history for a discretely valued character trait on the tree. The resultant tree is stored as a modified \code{"phylo"} object in stochastic character map (e.g., \code{\link{make.simmap}}) format.
}
\value{
	A modified phylogenetic tree of class \code{"phylo"} (or a modified \code{"multiPhylo"} object, for \code{nsim} > 1) with the following additional elements:
	\item{maps}{a list of named vectors containing the times spent in each state on each branch, in the order in which they occur.}
	\item{mapped.edge}{a matrix containing the total time spent in each state along each edge of the tree.}
	\item{states}{a vector containing the tip states.}
	\item{node.states}{a matrix containing the states at internal & terminal nodes (according to the dimensions of \code{edge}).}
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{plotSimmap}}, \code{\link{sim.rates}}
}
\keyword{phylogenetics}
\keyword{simulation}
