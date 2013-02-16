\name{splitplotTree}
\alias{splitplotTree}
\title{Plots a phylogeny in two columns}
\usage{
	splitplotTree(tree, fsize=1.0, ftype="reg", lwd=2, split=NULL, new.window=FALSE)
}
\arguments{
 	\item{tree}{an object of class "phylo".}
 	\item{fsize}{relative font size for tip labels.}
	\item{ftype}{font type - options are \code{"reg"}, \code{"i"} (italics), \code{"b"} (bold), or \code{"bi"} (bold-italics).}
	\item{lwd}{line width for plotting.}
	\item{split}{relative vertical position for splitting the tree (between 0 & 1).}
	\item{new.window}{whether or not to plot the split tree in a new window. If \code{FALSE} then the tree will be plotted in two columns within the same plotting window.}
}
\description{
	Function plots a tree in two columns or windows.
}
\value{
	plots a tree.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plotTree}}, \code{\link{plotSimmap}} 
}
\keyword{phylogenetics}
\keyword{plotting}
