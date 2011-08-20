\name{nodeHeights}
\alias{nodeHeights}
\title{Estimate diversity at each node of the tree}
\usage{
nodeHeights(tree)
}
\arguments{
  \item{tree}{is a phylogenetic tree in \code{"phylo"} format.}
}
\description{
  This function computes node heights for \code{\link{estDiversity}}.
}
\details{
}
\value{
  a matrix of the same dimensions as \code{tree$edge} containing the height above the root of each node in \code{edge}.
}
\references{
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{estDiversity}}, \code{\link{fitDiversityModel}}
}
\examples{
}
\keyword{phylogenetics}
