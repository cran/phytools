\name{add.species.to.genus}
\alias{add.species.to.genus}
\alias{genus.to.species.tree}
\title{Add species to genus on a phylogeny or bind simulated species subtrees to a backbone genus tree}
\usage{
add.species.to.genus(tree, species, genus=NULL, where=c("root","random"))
genus.to.species.tree(tree, species)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}. In the case of \code{genus.to.species.tree} this should be a genus-level backbone tree.}
	\item{species}{string contain species name in the format \code{"Genus_species"} or \code{"Genus species"}.}
	\item{genus}{for \code{add.species.to.genus}, optional argument containing the genus to which \code{species} is to be attached. If \code{NULL} then \code{genus} will be extracted from \code{species}.}
	\item{where}{for \code{add.species.to.genus}, the location to attach \code{species} to the tree. \code{where=} \code{"root"} will cause the species to be attached to the MRCA of all members of \code{genus}. \code{where="random"} will cause \code{species} to be attached at random to the subtree descended from the MRCA of all members of \code{genus}.}
}
\description{
	\code{add.species.to.genus} adds an additional species to a genus on a phylogeny. 
	
	\code{genus.to.species.tree} simulates pure-birth subtrees and then binds them at a random height along the terminal edge leading to each corresponding genus on a genus-level backbone tree.
}
\details{
	For \code{add.species.to.genus}, if \code{genus} contains only one species and \code{where="root"}, then \code{species} will be attached midway along the branch leading to the one species. If \code{where="random"} then \code{species} will be added at a random position along the edge. 
	
	If \code{genus} cannot be found in the tree, then the original tree is returned and a warning printed. If the tree is not ultrametric, then the resultant tree may not contain branch lengths and a warning will be printed. 
	
	Note that for some cases in which a tree is read in from file, R may initially think it is ultrametric, but then later (as tips are added) decide that it is not due to rounding of the edge lengths when it was written to file. This can most likely be resolved by using \code{\link{force.ultrametric}} to coerce the tree to be exactly ultrametric (to the maximum numerical precision permitted by R) before adding tips to the tree.
	
	If \code{genus} is non-monophyletic then \code{species} will be attached to the most inclusive group containing members of \code{genus} and a warning will be printed.
}
\value{
	An object of class \code{"phylo"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{add.random}}, \code{\link{bind.tip}}
}
\keyword{phylogenetics}
\keyword{utilities}
