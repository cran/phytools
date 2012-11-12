\name{ancThresh}
\alias{ancThresh}
\title{Ancestral character estimation under the threshold model using Bayesian MCMC}
\usage{
ancThresh(tree, x, ngen=1000, sequence=NULL, method="mcmc", control=list())
}
\arguments{
	\item{tree}{phylogenetic tree.}
	\item{x}{named vector containing discrete character states.}
	\item{ngen}{number of generations to run the MCMC.}
	\item{sequence}{assumed ordering of the discrete character state. If not supplied then numerical/alphabetical order is assumed!}
	\item{method}{only method currently available is \code{"mcmc"}.}
	\item{control}{list containing the following elements: \code{sample}, the sampling interval; \code{propliab}variance of the proposal distribution for liabilities; \code{propthresh} variance on the proposal distribution for the thresholds; \code{pr.anc} prior probability distribution on the ancestral states for each node, in a matrix - not all nodes need to be supplied; \code{pr.th} prior density on the thresholds; \code{burnin} number of generations to exclude for burn-in when plotting posterior probabilities on the tree; \code{plot} logical value indicating whether or not to plot the posterior probabilities; \code{print} logical value indicating whether or not to print the state of the MCMC; and \code{piecol} colors for the posterior probabilities plotted as pie charts at internal nodes.}
}
\description{
	This function uses Bayesian MCMC to estimate ancestral states and thresholds for a discrete character under the threshold model from quantitative genetics (Felsenstein 2012).
}
\value{
	This function returns as list with four elements:
	\item{ace}{posterior probabilities for each character state at each internal node.}
	\item{mcmc}{full posterior sample for the states.}
	\item{par}{full posterior sample for the thresholds and likelihood.}
	\item{liab}{full posterior sample of the liabilities at internal and tip nodes.}
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{anc.Bayes}}, \code{\link{threshBayes}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{bayesian}
