\name{plotTree}
\alias{plotTree}
\title{Plots rooted phylogenetic tree}
\usage{
plotTree(tree,color=NULL,fsize=1.0,ftype="reg",lwd=2,pts=TRUE)
}
\arguments{
  \item{tree}{a phylogenetic tree in \code{"phylo"} format; or multiple trees as an object of class \code{"multiPhylo"}.}
  \item{color}{color for the branches of the tree.}
  \item{fsize}{a numeric value specifying relative font size for tip labels.}
  \item{ftype}{a string specifying which font type to use. Options are: \code{"reg"},\code{"b"},\code{"i"}, \code{"bi"}, and \code{"off"}.}
  \item{lwd}{a numeric value specifying edge width.}
  \item{pts}{a logical value indicating whether to plot points at nodes.}
}
\description{
  This function plots a rooted phylogram.
}
\details{
}
\value{
  This function plots a rooted phylogram.
}
\references{
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{plot.phylo}}, \code{\link{plotSimmap}}
}
\examples{
tree<-pbtree(n=25)
plotTree(tree,color="blue",ftype="i",pts=FALSE)
}
\keyword{phylogenetics}
\keyword{plotting}
