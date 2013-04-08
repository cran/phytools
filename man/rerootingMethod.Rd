\name{rerootingMethod}
\alias{rerootingMethod}
\title{Get marginal ancestral state reconstructions by re-rooting}
\usage{
	rerootingMethod(tree, x, model=c("ER","SYM"))
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a vector of tip values for species; \code{names(x)} should be the species names.}
	\item{model}{any revsersible model. \code{model=c("ER","SYM")} recommended.}
}
\description{
	This function uses the re-rooting method of Yang et al. (1995) to get the marginal ancestral state estimates for each internal node of the tree using likelihood. This method get the conditional scaled likelihoods for the root node (which is the same as the marginal ancestral state reconstruction for that node) and successively moves the root to each node in the tree.
}
\value{
	a list containing the following elements:
	\item{loglik}{the log-likelihood.}
	\item{Q}{the fitted transition matrix between states.}
	\item{marginal.anc}{the marginal ancestral state reconstructions for each node.}
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ace}}, \code{\link{make.simmap}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
