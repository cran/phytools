\name{map.to.singleton}
\alias{map.to.singleton}
\alias{plotTree.singletons}
\title{Converts a mapped tree to a tree with singleton nodes}
\usage{
map.to.singleton(tree)
plotTree.singletons(tree)
}
\arguments{
	\item{tree}{a modified object of class \code{"phylo"} with a mapped discrete character or (for \code{plotTree.singletons} a tree with singleton nodes.}
}
\description{
	The function \code{map.to.singleton} takes an object of class \code{"phylo"} with a mapped discrete character and converts it to a tree with singleton nodes, in which edge has only one state. The states for each edge are stored in \code{names(tree$edge.length)}. \code{plotTree.singletons} plots a tree with singleton nodes.
}
\value{
	An object of class \code{"phylo"} with singleton nodes. \code{plotTree.singletons} plots a tree. If \code{names(tree$edge.length)!=NULL} it will use a different color from \code{\link{palette}} for each mapped state.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{collapse.singles}}, \code{\link{make.simmap}}
}
\keyword{phylogenetics}
\keyword{utilities}
