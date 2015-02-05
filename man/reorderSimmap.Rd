\name{reorderSimmap}
\alias{reorderSimmap}
\title{Reorder edges of a simmap tree}
\usage{
reorderSimmap(tree, order="cladewise")
}
\arguments{
	\item{tree}{a modified object of class \code{"phylo"}.}
	\item{order}{either \code{"cladewise"} or \code{"pruningwise"} (see \code{\link{reorder.phylo}}).}
}
\description{
	Function returns a reordered modified \code{"phylo"} object by using \code{reorder.phylo} but then sorting the additional elements \code{$mapped.edge} and \code{$maps} to have the same order as \code{$edge}.
}
\value{
	A modified object of class \code{"phylo"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{reorder.phylo}}, \code{\link{plotSimmap}}
}
\keyword{phylogenetics}
\keyword{utilities}
\keyword{comparative method}
