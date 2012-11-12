\name{fastMRCA}
\alias{fastMRCA}
\title{Get the MRCA of a pair of tip taxa}
\usage{
fastMRCA(tree, sp1, sp2)
}
\arguments{
  \item{tree}{a phylogenetic tree as an object of class \code{"phylo"}.}
  \item{sp1}{species name.}
  \item{sp2}{species name.}
}
\description{
	This function returns the most recent common ancestor (node number) for a pair of taxa. 
}
\details{
	This function is mostly redundant with \code{\link{findMRCA}} but for very large trees will be considerably faster. Function uses \code{\link{Ancestors}} in the phangorn package.
}
\value{
	The node number of the MRCA.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{Ancestors}}, \code{\link{findMRCA}}, \code{\link{mrca}}
}
\examples{
	tree<-pbtree(n=2000)
	anc<-fastMRCA(tree,"t1","t15")
}
\keyword{phylogenetics}
\keyword{utilities}
