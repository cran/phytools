\name{plotTree}
\alias{plotTree}
\title{Plots rooted phylogenetic tree}
\usage{
plotTree(tree, ...)
}
\arguments{
  \item{tree}{a phylogenetic tree in \code{"phylo"} format; or multiple trees as an object of class \code{"multiPhylo"}.}
  \item{...}{optional arguments.}
}
\description{
  This function plots a rooted phylogram. Arguments in \code{...} are passed to \code{\link{plotSimmap}}, with the exception of optional argument \code{color} which is used to determine the plotted color of the branch lengths of the tree. Note that the optional argument \code{pts}, which is by default \code{pts=TRUE} in \code{plotSimmap} has the default value of \code{pts=FALSE} in \code{plotTree}.
}
\value{
  This function plots a rooted phylogram.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{plotSimmap}}
}
\examples{
tree<-pbtree(n=25)
plotTree(tree,color="blue",ftype="i")
}
\keyword{phylogenetics}
\keyword{plotting}
