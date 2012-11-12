\name{phylomorphospace3d}
\alias{phylomorphospace3d}
\title{Creates phylomorphospace plot}
\usage{
phylomorphospace3d(tree, X, A=NULL, label=TRUE, control=list())
}
\arguments{
  \item{tree}{a phylogenetic tree in \code{"phylo"} format.}
  \item{X}{an \code{n} x 3 matrix of tip values for two characters in \code{n} species.}
  \item{A}{an optional \code{m} x 3 matrix (for \code{m} nodes) of values for two taits at internal nodes in the tree - if not supplied, these values will be estimated using \code{\link{anc.ML}}.}
  \item{label}{logical value indicating whether to print tip labels next to terminal nodes in the plot (presently doesn't do anything, but labels can be dropped using \code{control}).}
spin=TRUE,axes=TRUE,box=TRUE,simple.axes=FALSE,lwd=1,ftype="reg"
  \item{control}{a list containing the following optional control parameters: \code{spin}: a logical value indicating whether to animate the plot when created; \code{axes}: a logical indicating whether to plot the axes; \code{box}: a logical value indicating whether to plot in box; \code{simple.axes}: logical value indicating whether to replace \code{box} and \code{axes} with simpler axes; \code{lwd}: line widths; \code{ftype}: font type ("off" turns off labels altogether).}
}
\description{
  This function creates a phylomorphospace plot for three characters using the 3D visualization package, 'rgl'.
}
\value{
  this function creates a three dimensional phylomorphospace plot
}
\references{
   Sidlauskas, B. 2008. Continuous and arrested morphological diversification in sister clades of characiform fishes: A phylomorphospace approach. \emph{Evolution}, \bold{62}, 3135--3156.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
   tree<-pbtree(n=25)
   X<-fastBM(tree,nsim=3)
   phylomorphospace3d(tree,X,control=list(spin=FALSE))
}
\seealso{
  \code{\link{fancyTree}}, \code{\link{phenogram}}, \code{\link{phylomorphospace}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{plotting}
