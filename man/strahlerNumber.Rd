\name{strahlerNumber}
\alias{strahlerNumber}
\alias{extract.strahlerNumber}
\title{Computes Strahler number for trees and nodes}
\usage{
strahlerNumber(tree, plot=TRUE)
extract.strahlerNumber(tree, i, plot=TRUE)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{i}{order of Strahler number to extract for \code{extract.strahlerNumber}.}
	\item{plot}{logical value indicating whether to plot the tree with Strahler numbers for node labels.}
}
\description{
	Computes the Strahler number of all nodes and tips in a phylogenetic tree.
}
\details{
	The function \code{strahlerNumber} computes the Strahler number of all nodes and tips in the tree. For more information about Strahler numbers see \url{https://en.wikipedia.org/wiki/Strahler_number}. The function \code{extract.strahlerNumber} extracts all of the most inclusive clades of Strahler number \code{i}.
}
\value{
	Either a vector with the Strahler number for each tip and internal node; or (for \code{extract.strahlerNumber} the set of (most inclusive) subtrees with Strahler number \code{i} as an object of class \code{"multiPhylo"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
