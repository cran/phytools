\name{ltt}
\alias{ltt}
\title{Creates lineage-through-time plot (including extinct lineages)}
\usage{
ltt(tree,plot=TRUE,drop.extinct=FALSE,log.lineages=TRUE,gamma=TRUE)
}
\arguments{
  \item{tree}{is a phylogenetic tree in \code{"phylo"} format.}
  \item{plot}{a logical value indicating whether or not to create LTT plot.}
  \item{drop.extinct}{logical value indicating whether or not to drop extinct tips from the tree.}
  \item{log.lineages}{logical value indicating whether LTT plot should be on log-linear (default) or linear-linear scale.}
  \item{gamma}{logical value indicating whether or not to compute gamma from Pybus & Harvey (2000; \emph{Proc. Roy. Soc. B}).}
}
\description{
  This function computes LTT plot with extant and extinct lineages, and optionally conducts gamma-test of Pybus & Harvey (2000).
}
\details{
}
\value{
  a list with the following components:
  \item{times}{a vector of branching times.}
  \item{ltt}{a vector of linages.}
  \item{gamma}{optionally, a value for the gamma-statistic.}
  \item{p}{two-tailed P-value for the gamma-test.}
}
\references{
  Pybus, O. G., and P. H. Harvey. 2000. Testing macro-evolutionary models using incomplete molecular phylogenies. \emph{Proc. R. Soc. Lond. B}, \bold{267}, 2267-2272.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{gammatest}}
}
\examples{
}
\keyword{phylogenetics}
