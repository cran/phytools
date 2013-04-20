\name{make.simmap}
\alias{make.simmap}
\title{Simulate stochastic character maps on a phylogenetic tree or trees}
\usage{
	make.simmap(tree, x, model="SYM", nsim=1, ...)
}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"}, or a list of trees as an object of class \code{"multiPhylo"}.}
	\item{x}{a vector containing the tip states for a discretely valued character, or a matrix containing the prior probabilities of tip states in rows.}
	\item{model}{a character string containing the model - options as in \code{\link{ace}}.}
	\item{nsim}{number of simulations. If \code{tree} is an object of class \code{"multiPhylo"}, then \code{nsim} simulations will be conducted per tree.}
	\item{...}{optional arguments. So far, \code{pi} gives the state frequencies - options are \code{"equal"}, \code{"estimated"}, or a vector with the frequencies. If \code{pi="estimated"} then the stationary distribution is estimated by numerically solving \code{pi*Q=0} for \code{pi}. Defaults to \code{pi="equal"} which results in the root node being sampled from the conditional scaled likelihood distribution at the root. \code{message} tells whether or not to print a message containing the rate matrix, \emph{Q} and state frequencies. Defaults to \code{message=TRUE}. \code{tol} gives the tolerance for zero elements in \code{Q}. (Elements less then \code{tol} will be reset to \code{tol}.}
}
\description{
	This function fits a continuous-time reversible Markov model for the evolution of \code{x} and then simulates stochastic character histories using that model and the tip states on the tree. This is the same procedure that is described in Bollback (2006), except that simulation is performed using a fixed value of \emph{Q} instead of by sampling \emph{Q} from a user-specified prior distribution. 
}
\details{
	Samples trees from the conditional posterior distribution - conditioning on the maximum likelihood transition matrix given our model. Uses code modified from \code{\link{ace}} by Paradis et al.

	As of phytools>=0.2-33 \code{x} can be a vector of states or a matrix containing the prior probabilities of tip states in rows. In this case the column names of \code{x} should contain the states, and the row names should contain the tip names.

	Note that there was a small (but potentially significant) bug in how node states were simulated by \code{make.simmap} in versions of phytools<=0.2-26. Between phytools 0.2-26 and 0.2-36 there was also a bug for asymmetric models of character change (e.g., \code{model="ARD"}). Finally, between phytools 0.2-33 and phytools 0.2-47 there was an error in use of the conditional likelihoods for the root node, which caused the root node of the tree to be sampled incorrectly. All of these issues should be fixed in the present version.
}
\value{
 	A modified phylogenetic tree of class \code{"phylo"} (or a modified \code{"multiPhylo"} object, for \code{nsim} > 1) with the following additional elements:
 	\item{maps}{a list of named vectors containing the times spent in each state on each branch, in the order in which they occur.}
 	\item{mapped.edge}{a matrix containing the total time spent in each state along each edge of the tree.}
}
\references{
 	Paradis, E., J. Claude, and K. Strimmer. 2004. APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289--290.

	Huelsenbeck, J. P., R. Neilsen, and J. P. Bollback. 2003. Stochastic mapping of morphological characters. \emph{Systematic Biology}, \bold{52}, 131-138. 

	Bollback, J. P. 2006. Stochastic character mapping of discrete traits on phylogenies. \emph{BMC Bioinformatics}, \bold{7}, 88.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownie.lite}}, \code{\link{brownieREML}}, \code{\link{evol.vcv}}, \code{\link{read.simmap}}, \code{\link{plotSimmap}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{simulation}
\keyword{bayesian}

