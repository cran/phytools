\name{add.color.bar}
\alias{add.color.bar}
\title{Add color bar to a plot}
\usage{
add.color.bar(leg, cols, title=NULL, lims=c(0,1), digits=1, prompt=TRUE, ...)
}
\arguments{
	\item{leg}{numerical value for the length of the legend.}
	\item{cols}{colors for the legend.}
	\item{title}{text to plot above the bar.}
	\item{lims}{range for the bar.}
	\item{digits}{digits for plotted numbers.}
	\item{prompt}{logical value indicating whether the location of the legend should be obtained interactively.}
	\item{...}{optional arguments including: \code{x} x-coordinate of the legend (if \code{prompt=FALSE}); \code{y} y-coordinate of the legend; \code{subtitle} optional legend subtitle.}
}
\description{
	This function adds a color bar to a plot created by \code{\link{plotBranchbyTrait}}. The only way to add the color bar at present is by clicking on the desired location within your plot. It is also used internally by \code{contMap(...,type="fan")} and \code{densityMap(...,type="fan")}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{contMap}}, \code{\link{densityMap}}, \code{\link{plotBranchbyTrait}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
