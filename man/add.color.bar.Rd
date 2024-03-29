\name{add.color.bar}
\alias{add.color.bar}
\title{Add color bar to a plot}
\usage{
add.color.bar(leg, cols, title=NULL, lims=c(0,1), digits=1, prompt=TRUE,
   lwd=4, outline=TRUE, ...)
}
\arguments{
	\item{leg}{numerical value for the length of the legend.}
	\item{cols}{colors for the legend.}
	\item{title}{text to plot above the bar.}
	\item{lims}{range for the bar.}
	\item{digits}{digits for plotted numbers.}
	\item{prompt}{logical value indicating whether the location of the legend should be obtained interactively.}
	\item{lwd}{width of the plotted bar.}
	\item{outline}{logical value indicated whether or not to outline the plotted color bar with a 1 pt line.}
	\item{...}{optional arguments including: \code{x} x-coordinate of the legend (if \code{prompt=FALSE}); \code{y} y-coordinate of the legend; \code{subtitle} optional legend subtitle; \code{direction} direction of the color bar (i.e., increase from left to right or from right to left); and \code{fsize}, font size for the legend text.}
}
\description{
	This function adds a color gradient legend to an existing plot.
}
\details{
	This function was designed to add a color bar to a plot created by \code{\link{plotBranchbyTrait}}. A color bar can be added by clicking on a location within the plot (when \code{prompt=TRUE}) or by setting \code{prompt=FALSE} and supplying x & y coordinates for the object. 
	
	The function is also used internally by S3 methods \code{\link{plot.contMap}} and \code{\link{plot.densityMap}}, as well as by \code{\link{errorbar.contMap}}.
}
\references{
	Revell, L. J. (2013) Two new graphical methods for mapping trait evolution on phylogenies. \emph{Methods in Ecology and Evolution}, \bold{4}, 754-759.
	
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{contMap}}, \code{\link{densityMap}}, \code{\link{errorbar.contMap}}, \code{\link{plotBranchbyTrait}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
\keyword{utilities}
