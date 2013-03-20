\name{repPhylo}
\alias{repPhylo}
\title{Replicate a tree into a list of trees}
\usage{
	repPhylo(tree, times)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{times}{number of times to replicate tree.}
}
\description{
	This function is functionally equivalent to \code{\link{rep}}, except that it works for phylogenies. Will create an object of class \code{"multiPhylo"}.
}
\value{
	An object of class \code{"multiPhylo"}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
	tree<-pbtree(n=100)
	trees<-repPhylo(tree,100)
}
\seealso{
	\code{\link{c.phylo}}, \code{\link{rep}}
}
\keyword{phylogenetics}
\keyword{utilities}
