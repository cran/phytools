\name{xkcdTree}
\alias{xkcdTree}
\title{Plot xkcd style tree}
\usage{
xkcdTree(tree, file=NULL, gsPath=NULL, fsize=2, lwd=7, color="blue",
     dim=c(8.5,11), jitter=0.001, waver=c(0.1,0.1), tilt=0, right=TRUE)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{file}{an output file name. If the ending ".pdf" is not provided, it will be automatically appended.}
	\item{gsPath}{for Windows, the full path to the Ghostscript executable. If not provided, the default for gs9.06 is used.}
	\item{fsize}{font size.}
	\item{lwd}{line width.}
	\item{color}{color for lines.}
	\item{dim}{vector of horizontal and vertical dimensions for the output PDF.}
	\item{jitter}{amount of stochasticity ("jitter") for plotting - higher is more jittery.}
	\item{waver}{a vector with approximate horizontal and vertical wavelength for plotting - lower is more wavery.}
	\item{tilt}{declination from horizontal (or vertical, if \code{right=FALSE} for tip labels, in degrees.}
	\item{right}{logical indicating whether to plot the tree with rightward orientation (if \code{TRUE}) or upwards.}
}
\description{
	Function plots a tree in the style of \url{http://xkcd.com} (for instance, \url{http://xkcd.com/867/}).
}
\details{
	Uses package "extrafont" and font "xkcd.ttf". Before using, install package "extrafont" and add font "xkcd.ttf". Hopefully, if you have not done that, a sensible error will still be returned.
}
\value{
	plots a tree to PDF.
}
\author{
	Liam Revell \email{liam.revell@umb.edu}; an anonymous user also contributed much of the code for internally called function \code{xkcdLine}.
}
\seealso{
	\code{\link{fancyTree}}, \code{\link{plot.phylo}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
