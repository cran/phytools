\name{phyl.vcv}
\alias{phyl.vcv}
\title{Compute evolutionary VCV matrix for a tree & dataset}
\usage{
phyl.vcv(X,C,lambda)
}
\arguments{
  \item{lambda}{value for lambda transformation.}
  \item{X}{data matrix.}
  \item{C}{matrix containing the height above the root of each pair of species in the tree. Typically this will have been produced by calling \code{\{link{vcv.phylo}}.}
}
\description{
  Primarily an internal function for \code{\link{phyl.pca}}; this can be used to compute the phylogenetic trait variance-covariance matrix given a phylogenetic VCV, lambda, and a data matrix.
}
\details{
  Do not use unless you know what you're doing.
}
\value{
  matrix.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{statistics}
\keyword{utilities}
