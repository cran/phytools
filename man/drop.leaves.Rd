\name{drop.leaves}
\alias{drop.leaves}
\title{Drop all the leaves (tips) from a tree}
\usage{
drop.leaves(tree, ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{...}{optional arguments. Presently includes only the logical value \code{keep.tip.labels} which tells the function how to labels the tips on the reduced tree.}
}
\description{
	Drops all the leaves from a tree, leaving behind only the structure leading to internal nodes.
}
\value{
	An object of class \code{"phylo"}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
