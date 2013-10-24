\name{plot.backbonePhylo}
\alias{plot.backbonePhylo}
\title{Plots backbone tree with triangles as clades}
\usage{
\method{plot}{backbonePhylo}(x, ...)
}
\arguments{
	\item{x}{an object of class \code{"backbonePhylo"}.}
	\item{...}{optional arguments. So far includes only \code{vscale}, which is used to rescale the vertical dimension in plotting.}
}
\description{
	Function plots a backbone tree (stored as an object of class \code{"backbonePhylo"}) with triangles as subtrees.
}
\value{
	Plots a tree.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{phylo.toBackbone}}
}
\keyword{phylogenetics}
\keyword{plotting}

