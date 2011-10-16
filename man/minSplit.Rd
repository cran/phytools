\name{minSplit}
\alias{minSplit}
\title{Finding the minimum (median) split in the posterior sample}
\usage{
minSplit(tree,split.list,method="sum",printD=FALSE)
}
\arguments{
  \item{tree}{a phylogenetic tree in \code{"phylo"} format.}
  \item{split.list}{either a matrix with two named columns, \code{"node"} and \code{"bp"}; a \code{$mcmc} matrix from \code{evol.rate.mcmc()}; or the entire raw output from \code{evol.rate.mcmc()}.}
  \item{method}{an optional string indicating the criterion to minimize: options are \code{"sum"} and \code{"sumsq"}.}
  \item{printD}{logical specifying whether to print distances to screen (\code{FALSE} by default).}
}
\description{
  This function takes a phylogenetic tree and a list of splits and identifies the split with the smallest summed or summed squared distances to all the other splits.  Used to be called \code{min.split()}.
}
\details{
}
\value{
  a list with the following components:
  \item{node}{node for the minimum split.}
  \item{bp}{location on the branch leading to \code{node} of the minimum split.}
}
\references{
  Revell, L. J., D. L. Mahler, P. Peres-Neto, and B. D. Redelings. In review.
  A new method for identifying exceptional phenotypic diversification.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{evol.rate.mcmc}}, \code{\link{posterior.evolrate}}
}
\examples{
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{bayesian}
