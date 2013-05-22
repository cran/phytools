\name{nodeHeights}
\alias{nodeHeights}
\title{Compute the heights above the root of each node}
\usage{
nodeHeights(tree)
}
\arguments{
	\item{tree}{a phylogeny as an object of class \code{"phylo"}.}
}
\description{
	This function computes node heights.
}
\details{
	The function also gives a handy way to get the total length of the tree from the root to the heighest tip which will be given by \code{max(nodeHeights(tree))}.
}
\value{
	A matrix of the same dimensions as \code{tree$edge} containing the height above the root of each node in \code{edge}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{estDiversity}}, \code{\link{fitDiversityModel}}
}
\examples{
tree<-rtree(10)
X<-nodeHeights(tree)
}
\keyword{phylogenetics}
\keyword{utilities}
\keyword{comparative method}
