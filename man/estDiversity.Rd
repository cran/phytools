\name{estDiversity}
\alias{estDiversity}
\title{Estimate diversity at each node of the tree}
\usage{
estDiversity(tree,x)
}
\arguments{
  \item{tree}{is a phylogenetic tree in \code{"phylo"} format.}
  \item{x}{a vector containing the biogeographic area for each of the tip taxa.}
}
\description{
  This function estimates the lineage density at each node in the tree based on a biogeographic model (based on Mahler et al. 2010; \emph{Evolution}).
}
\details{
  Very slow. This function (presently) relies on many calls to \code{ace}.
}
\value{
  a vector containing the estimated lineage density at each node
}
\references{
  Mahler, D. L, L. J. Revell, R. E. Glor, and J. B. Losos. 2010. Ecological opportunity and the rate of morphological evolution in the diversification of Greater Antillean anoles. \emph{Evolution}, \bold{64}, 2731-2745.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{fitDiversityModel}}
}
\examples{
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
