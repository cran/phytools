\name{plotSimmap}
\alias{plotSimmap}
\title{Plot stochastic character mapped tree}
\usage{
	plotSimmap(tree, colors=NULL, fsize=1.0, ftype="reg", lwd=2, pts=TRUE, 
     node.numbers=FALSE, mar=NULL, add=FALSE, offset=NULL, direction="rightwards")
}
\arguments{
 	\item{tree}{a modified object of class "phylo" or "multiPhylo" containing a stochastic mapping or set of mappings (e.g., see \code{\link{read.simmap}} & \code{\link{make.simmap}}).}
 	\item{colors}{a vector with names translating the mapped states to colors - see \code{Examples}.}
 	\item{fsize}{relative font size for tip labels.}
	\item{ftype}{font type - options are \code{"reg"}, \code{"i"} (italics), \code{"b"} (bold), or \code{"bi"} (bold-italics).}
	\item{lwd}{line width for plotting.}
	\item{pts}{logical value indicating whether or not to plot filled circles at each vertex of the tree, as well as at transition points between mapped states.  Default is \code{TRUE}.}
	\item{node.numbers}{a logical value indicating whether or not node numbers should be plotted.}
	\item{mar}{vector containing the margins for the plot to be passed to \code{\link{par}}. If not specified, the default margins are [0.1,0.1,0.1,0.1].}
	\item{add}{a logical value indicating whether or not to add the plotted tree to the current plot (\code{TRUE}) or create a new plot (\code{FALSE}, the default).}
	\item{offset}{offset for the tip labels. Primarily to be used internally by \code{\link{densityMap}}.}
	\item{direction}{plotting direction. Options are \code{"rightwards"} (the default) and \code{"leftwards"}. Note that for some reason that is not totally clear, \code{node.numbers=TRUE} does not work for \code{direction="leftwards"}.}
}
\description{
	Function plots one or multiple stochastic character mapped trees.
}
\details{
	The underscore character \code{"_"} is automatically swapped for a space in tip labels, as in \code{\link{plot.phylo}}.
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
	\code{\link{densityMap}}, \code{\link{make.simmap}}, \code{\link{read.simmap}} 
}
\examples{
	Q<-matrix(c(-2,1,1,1,-2,1,1,1,-2),3,3)
	rownames(Q)<-colnames(Q)<-letters[1:3]
	tree<-sim.history(pbtree(n=100,scale=1),Q)
	cols<-setNames(c("blue","red","green"),letters[1:3])
	plotSimmap(tree,cols,ftype="i",fsize=0.7)
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
