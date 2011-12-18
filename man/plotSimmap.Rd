\name{plotSimmap}
\alias{plotSimmap}
\title{Plot stochastic character mapped tree}
\usage{
plotSimmap(tree,colors=NULL,fsize=1.0,ftype="reg",lwd=2,pts=TRUE)
}
\arguments{
  \item{tree}{a modified object of class "phylo" or "multiPhylo" containing a stochastic mapping or set of mappings (e.g., see \code{\link{read.simmap}} & \code{\link{make.simmap}}).}
  \item{colors}{a vector with names translating the mapped states to colors - see \code{Examples}.}
  \item{fsize}{relative font size for tip labels.}
  \item{ftype}{font type - options are \code{"reg"}, \code{"i"} (italics), \code{"b"} (bold), or \code{"bi"} (bold-italics).}
	\item{lwd}{line width for plotting.}
  \item{pts}{logical value indicating whether or not to plot filled circles at each vertex of the tree, as well as at transition points between mapped states.  Default is \code{TRUE}.}
}
\description{
  Function plots one or multiple stochastic character mapped trees.
}
\details{
}
\value{
  plots a tree.
}
\references{
  Huelsenbeck, J. P., R. Neilsen, and J. P. Bollback. 2003. Stochastic mapping of morphological characters. \emph{Systematic Biology}, \bold{52}, 131-138. 

  Bollback, J. P. 2006. Stochastic character mapping of discrete traits on phylogenies. \emph{BMC Bioinformatics}, \bold{7}, 88.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{read.simmap}}, \code{\link{make.simmap}}
}
\examples{
X<-"((A:{aqua,0.3:terr,0.4:aqua,0.3},B:{aqua,1.0}):{aqua,1.0},(C:{aqua,0.25:terr,0.75},D:{aqua,1.0}):{aqua,1.0});"
phy<-read.simmap(text=X)
leg<-c("blue","red"); names(leg)<-c("aqua","terr")
plotSimmap(phy,leg,ftype="i")
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
