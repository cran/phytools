\name{plotTree.wBars}
\alias{plotTree.wBars}
\title{Plot a tree with bars at the tips}
\usage{
plotTree.wBars(tree, x, scale=1, width=NULL, type="phylogram", method="plotTree", 
	...)
}
\arguments{
 	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a named vector of trait values (normally > 0).}
 	\item{scale}{scaling factor for the tip bars (relative to the total tree height).}
	\item{width}{width of the tip bars.}
	\item{type}{plot type. Can be \code{"phylogram"} or \code{"fan"}.}
	\item{method}{plotting method to use. Can be \code{"plotTree"} (for \code{\link{plotTree}}) or \code{"plotSimmap"} (for \code{\link{plotSimmap}}).}
	\item{...}{optional arguments to be passed to code{\link{plotTree}} or \code{\link{plotSimmap}}.}
}
\description{
	Function plots a phylogeny in phylogram or fan style with bars at the tips representing the values for a phenotypic trait.
}
\value{
	Plots a tree.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plotSimmap}}, \code{\link{plotTree}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}