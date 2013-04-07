\name{findMRCA}
\alias{findMRCA}
\title{Get the MRCA of a set of taxa}
\usage{
findMRCA(tree, tips=NULL)
}
\arguments{
  \item{tree}{a phylogenetic tree as an object of class \code{"phylo"}.}
  \item{tips}{a vector containing a set of tip labels.}
}
\description{
	This function returns the most recent common ancestor (node number) for a set of taxa.
}
\details{
	If \code{tips==NULL} will return the result of a normal function call to \code{\link{mrca}}.
}
\value{
  The node number of the MRCA, or a matrix of node numbers (if \code{tips==NULL}).
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{mrca}}
}
\examples{
	tree<-pbtree(n=20)
	anc<-findMRCA(tree,c("t1","t10","t15"))
}
\keyword{phylogenetics}
\keyword{utilities}
