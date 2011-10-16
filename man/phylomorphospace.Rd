\name{phylomorphospace}
\alias{phylomorphospace}
\title{Creates phylomorphospace plot}
\usage{
phylomorphospace(tree,X,A=NULL,label=TRUE,control=list())
}
\arguments{
  \item{tree}{a phylogenetic tree in \code{"phylo"} format.}
  \item{X}{an \code{n} x 2 matrix of tip values for two characters in \code{n} species.}
  \item{A}{an optional \code{m} x 2 matrix (for \code{m} nodes) of values for two taits at internal nodes in the tree - if not supplied, these values will be estimated using \code{\link{ace}}.}
  \item{label}{logical value indicating whether to print tip labels next to terminal nodes in the plot.}
  \item{control}{a list containing the following optional control parameters: \code{col.edge}: a vector of edge colors; and \code{col.node}: a vector of node colors.}
}
\description{
  This function creates a phylomorphospace plot for two characters following Sidlauskas (2008; \emph{Evolution}).
}
\details{
}
\value{
  this function creates a phylomorphospace plot
}
\references{
  Paradis, E., J. Claude, and K. Strimmer. 2004. APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289--290.

  Sidlauskas, B. 2008. Continuous and arrested morphological diversification in sister clades of characiform fishes: A phylomorphospace approach. \emph{Evolution}, \bold{62}, 3135--3156.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
}
\examples{
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{plotting}
