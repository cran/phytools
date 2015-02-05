\name{setMap}
\alias{setMap}
\title{Set color map for objects of class \code{"contMap"} or \code{"densityMap"}}
\usage{
setMap(x, ...)
}
\arguments{
	\item{x}{an object of class \code{"contMap"} or \code{"densityMap"}.}
	\item{...}{arguments to be passed to \code{\link{colorRampPalette}}.}
}
\description{
	Function to change the color map (ramp) in an object of class \code{"contMap"} or \code{"densityMap"}.
}
\value{
	An object of class \code{"contMap"} or \code{"densityMap"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Revell, L. J. 2013. Two new graphical methods for mapping trait evolution on phylogenies. \emph{Methods in Ecology and Evolution}, \bold{4}, 754-759.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{contMap}}, \code{\link{densityMap}}
}
\keyword{phylogenetics}
\keyword{plotting}

