\name{phenogram}
\alias{phenogram}
\title{Plot phenogram (traitgram)}
\usage{
	phenogram(tree, x, fsize=1.0, ftype="reg", colors=NULL, axes=list(),
		add=FALSE, ...)
}
\arguments{
	\item{tree}{an object of class "phylo", with or without a mapped discrete character.}
	\item{x}{a vector containing the states at the tips \emph{or} the states at all the tips and the internal nodes of the tree.}
	\item{fsize}{relative font size for tip labels.}
	\item{ftype}{font type - options are \code{"reg"}, \code{"i"} (italics), \code{"b"} (bold), or \code{"bi"} (bold-italics).}
	\item{colors}{colors for plotting the mapped character (if available) in \code{tree}. If no character is mapped on the tree, then a single color for all the branches of the tree can be provided.}
	\item{axes}{list of axis dimensions. Items are \code{time} and \code{trait}.}
	\item{add}{optional logical value indicating whether to add to an open plot. If \code{TRUE}, then new axes will not be plotted.}
	\item{...}{optional arguments including \code{xlim}, \code{ylim}, \code{log}, \code{main}, \code{sub}, \code{xlab}, \code{ylab}, \code{asp}, \code{type}, \code{lty}, \code{lwd}, \code{offset}, \code{digits}, and \code{digits} are as in \code{\link{plot.default}} or \code{\link{par}}. Note that \code{axes} overrides \code{xlim} and \code{ylim}. Finally, \code{spread.labels} is a logical value indicating whether or not to minimize tip label overlap (default is \code{FALSE}); \code{spread.cost} is a numeric vector indicating the relative penalty to be used for label overlap and deviance if \code{spread.labels=TRUE}; finally, \code{link} is a numeric value by which to offset the tip labels, linking them to the tips with a dashed line (default is \code{0} if \code{spread.labels=FALSE} or \code{0.1*max(nodeHeights(tree))} otherwise). The optional argument \code{offsetFudge} "fudges" the computation of label offset in scaling \code{xlim}. It is 1.37, which is the correct fudge in the Windows R GUI, but this may need to be changed in other systems.}
}
\description{
	Function plots a traitgram, that is, a projection of the phylogenetic tree in a space defined by phenotype (on the \emph{y} axis) and time (on the \emph{x}). If a discrete character is mapped on the tree this will also be plotted.
}
\details{
	For \code{spread.labels=TRUE} numerical optimization is performed to optimize the distribution of the labels vertically, where the solution depends on the vector \code{spread.cost} containing the cost of overlap (first) and the cost of deviation from the vertical position of the tip. Because this is done via numerical optimization, plotting may hang briefly while the best solution is found (especially for large trees).
}
\value{
	plots a traitgram, optionally with a mapped discrete character.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
	tree<-pbtree(n=20,scale=2)
	x<-fastBM(tree)
	phenogram(tree,x)
	# or, simulate a discrete character history
	tree<-sim.history(tree,Q=matrix(c(-1,1,1,-1),2,2),anc="1")
	# simulate in which the rate depends on the state
	x<-sim.rates(tree,c(1,10))
	phenogram(tree,x)
	# now use spread.labels
	tree<-pbtree(n=40)
	x<-fastBM(tree)
	phenogram(tree,x,spread.labels=TRUE,spread.cost=c(1,0))
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
