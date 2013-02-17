\name{roundBranches}
\alias{roundBranches}
\title{Rounds the branch lengths of a tree}
\usage{
	roundBranches(tree, digits)
}
\arguments{
	\item{tree}{tree.}
	\item{digits}{number of digits for rounding. Passed to \code{\link{round}}.}
}
\description{
	This function rounds the branch lengths of a tree and reconciles any mappings (as in \code{\link{read.simmap}}) with the rounded branch lengths.
}
\value{
  a tree with branch lengths, or modified \code{"phylo"} object with a mapped discrete character.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
