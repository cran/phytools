\name{multiC}
\alias{multiC}
\title{Returns a list with phylogenetic VCV matrix for each mapped state}
\usage{
multiC(tree)
}
\arguments{
	\item{tree}{a phylogeny with mapped discrete state in a modified object of class \code{"phylo"} (e.g., see \code{\link{read.simmap}}.}
}
\description{
	This function takes a modified \code{"phylo"} object as input and returns a set of so-called phylogenetic covariance matrices as a list: one for each mapped state.
}
\value{
	A list of matrices.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{evolvcv.lite}}, \code{\link{read.simmap}}, \code{\link{vcvPhylo}}, \code{\link{vcv.phylo}}
}
\keyword{phylogenetics}
\keyword{utilities}
