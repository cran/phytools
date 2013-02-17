\name{orderMappedEdge}
\alias{orderMappedEdge}
\title{Matches nodes between two trees}
\usage{
	orderMappedEdge(trees, ordering=NULL)
}
\arguments{
	\item{trees}{object of class \code{"phylo"} or \code{"multiPhylo"}.}
	\item{ordering}{ordering for the columns of \code{$mapped.edge}. If \code{NULL}, then an alphabetical order is assumed. Options are \code{"alphabetical"}, \code{"numerical"}, or any specific ordering of the mapped traits (e.g., \code{c("A","B","C")}.}
}
\description{
	This function takes a modified object of class \code{"multiPhylo"} with a mapped discrete character (e.g., see \code{\link{read.simmap}} and sorts the columns of each \code{tree$mapped.edge} to have the same state ordering. Also works for a single tree.
}
\value{
	a modified object of class \code{"phylo"} or \code{"multiPhylo"}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
