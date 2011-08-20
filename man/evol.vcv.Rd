\name{evol.vcv}
\alias{evol.vcv}
\title{Likelihood test for variation in the evolutionary VCV matrix}
\usage{
evol.vcv(scm.tre,dat,maxit=2000)
}
\arguments{
  \item{scm.tre}{a phylogenetic tree in modified \code{"phylo"} format (see \code{\link{read.simmap}}).}
  \item{dat}{an \code{n} x \code{m} matrix of tip values for \code{m} continuously valued traits in \code{n} species - row names should be species names.}
  \item{maxit}{an optional integer value indicating the maximum number of iterations for optimization - may need to be increased for large trees.}
}
\description{
  This function takes a modified \code{"phylo"} object with a mapped binary or multistate trait and data for an arbitrary number of continuously valued character.  It then fits the multiple evolutionary variance-covariance matrix (rate matrix) model of Revell & Collar (2009; \emph{Evolution}.
}
\details{
  May have to increase \code{maxit} for large trees and more than 2 traits.
}
\value{
  a list with the following components:
  \item{R.single}{vcv matrix for the single rate matrix model.}
  \item{logL1}{log-likelihood for single matrix model.}
  \item{k1}{number of parameters in the single marix model.}
  \item{R.multiple}{\code{m} x \code{m} x \code{p} array containing the \code{p} estimated vcv matrices for the \code{p} regimes painted on the tree.}
  \item{logL.multiple}{log-likelihood of the multi-matrix model.}
  \item{k2}{number of parameters estimated in this model.}
  \item{P.chisq}{P-value of the X^2 test on the likelihood ratio.}
  \item{convergence}{logical value indicating whether or not the optimization has converged.}
}
\references{
  Revell, L. J., and D. C. Collar. 2009. Phylogenetic analysis of the evolutionary correlation using likelihood. \emph{Evolution}, \bold{63}, 1090-1100.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{evol.rate.mcmc}}, \code{\link{brownie.lite}}
}
\examples{
}
\keyword{phylogenetics}
