\name{phyl.RMA}
\alias{phyl.RMA}
\title{Phylogenetic reduced major axis (RMA) regression}
\usage{
phyl.RMA(x,y,tree,method="BM")
}
\arguments{
  \item{x}{vector with names.}
  \item{y}{vector with names.}
  \item{tree}{a phylogenetic tree in \code{"phylo"} format.}
  \item{method}{method to obtain the correlation structure: can be \code{"BM"} or \code{"lambda"}.}
}
\description{
  This function performs phylogenetic RMA regression.
}
\details{
  Optionally jointly estimates lambda if \code{method="lambda"}. Likelihood optimization of lambda is performed using \code{\link{optimize}} on the interval (0,1).
}
\value{
  a list with the following components:
  \item{RMA.beta}{a vector of RMA regression coefficients.}
  \item{V}{a VCV matrix for the traits.}
  \item{lambda}{fitted value of lambda (\code{method="lambda"} only).}
  \item{logL}{log-likelihood (\code{method="lambda"} only).}
  \item{resid}{a vector of residuals for \code{y} given \code{x}.}
}
\references{
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{phyl.pca}}, \code{\link{phyl.resid}}
}
\examples{
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{statistics}
