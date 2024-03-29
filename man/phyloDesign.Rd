\name{phyloDesign}
\alias{phyloDesign}
\title{Compute design matrix for least squares analyses}
\usage{
phyloDesign(tree)
}
\arguments{
	\item{tree}{phylogenetic tree.}
}
\description{
	Primarily an internal function for \code{\link{optim.phylo.ls}}, this function creates a design matrix for least squares phylogenetic analysis.
}
\details{
    This function returns a matrix containing the edges in the tree (in columns) and pairs of tip node numbers (in rows). Values in the matrix are either \code{1} if the edge is on the shortest path between the two tips; and \code{0} otherwise. Probably do not use unless you know what you're doing.
}
\value{
	A matrix.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{phylogeny inference}
\keyword{least squares}
