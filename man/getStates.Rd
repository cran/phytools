\name{getStates}
\alias{getStates}
\title{Get the states at nodes or tips from a mapped tree}
\usage{
getStates(tree, type=c("nodes","tips","both"))
}
\arguments{
	\item{tree}{is a modified object of class \code{"phylo"} or \code{"multiPhylo"}.}
	\item{type}{mode indicating whether to get states at the nodes (\code{type="nodes"}) or the tips (\code{type="tips"}) of the tree.}
}
\description{
	Gets the states from the nodes or tips of a mapped tree (e.g., \code{\link{make.simmap}}).
}
\value{
	A named vector (for \code{"phylo"}) or matrix (for \code{"multiPhylo"}).
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{describe.simmap}}, \code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{sim.history}}
}
\keyword{phylogenetics}
\keyword{utilities}
\keyword{discrete character}
