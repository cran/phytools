\name{sim.rates}
\alias{sim.rates}
\title{Brownian simulationn with multiple evolutionary rates}
\usage{
sim.rates(mtree,sig2,anc=0,nsim=1,internal=F,plot=F)
}
\arguments{
  \item{mtree}{is a stochastic map format phylogenetic tree in modified \code{"phylo"} format (e.g., see \code{\link{make.simmap}}).}
  \item{sig2}{a named vector containing the rates for each state; names should be states in \code{mtree}.}
  \item{anc}{optional value for the root state.}
  \item{nsim}{number of simulations.}
  \item{internal}{logical value indicating whether to return states at internal nodes.}
  \item{plot}{logical value indicating whether or not to visual the rate heterogeneity (default value is \code{FALSE}.}
}
\description{
  This function conducts BM simulation on a tree with multiple rates.
}
\details{
}
\value{
  a vector (for \code{nsim=1}) or matrix containing the tip states for the \code{n} species in the tree.
}
\references{
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{fastBM}}, \code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{sim.history}}
}
\examples{
}
\keyword{phylogenetics}
\keyword{simulation}
\keyword{comparative method}
