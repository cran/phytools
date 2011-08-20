\name{phyl.pca}
\alias{phyl.pca}
\title{Phylogenetic principal components analysis}
\usage{
phyl.pca(tree,Y,method="BM",mode="cov")
}
\arguments{
  \item{tree}{a phylogenetic tree in \code{"phylo"} format.}
  \item{Y}{data matrix with traits in columns.}
  \item{method}{method to obtain the correlation structure: can be \code{"BM"} or \code{"lambda"}.}
  \item{mode}{is the mode for the PCA: can be \code{"cov"} or \code{"corr"}.}
}
\description{
  This function performs phylogenetic PCA (e.g., Revell 2009; \emph{Evolution}).
}
\details{
}
\value{
  a list with the following components:
  \item{Eval}{diagonal matrix of eigenvalues.}
  \item{Evec}{matrix with eigenvectors in columns.}
  \item{S}{matrix with scores.}
  \item{L}{matrix with loadings.}
  \item{lambda}{fitted value of lambda (\code{method="lambda"} only).}
  \item{logL}{log-likelihood for lambda model (\code{method="logL"} only).}
}
\references{
  Revell, L. J. 2009. Size-correction and principal components for interspecific comparative studies. \emph{Evolution}, \bold{63},
  3258--3268.

}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{phyl.resid}}
}
\examples{
}
\keyword{phylogenetics}
