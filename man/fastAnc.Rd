\name{fastAnc}
\alias{fastAnc}
\title{Fast estimation of ML ancestral states}
\usage{
fastAnc(tree, x)
}
\arguments{
  \item{tree}{an object of class \code{"phylo"}.}
  \item{x}{a vector of tip values for species; \code{names(x)} should be the species names.}
}
\description{
  This function performs fast estimation of the ML ancestral states for a continuous trait by taking advantage of the fact that the state computed for the root node of the tree during Felsenstein's (1985) contrasts algorithm is also the MLE of the root node. Thus, the function reroots the tree at all internal nodes and computes the contrasts state at the root each time.
}
\value{
  a named vector containing the states at internal nodes - names are node numbers.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{ace}}, \code{\link{anc.Bayes}}, \code{\link{anc.ML}}, \code{\link{pic}}
}
\examples{
tree<-rtree(50)
x<-fastBM(tree) # simulate using fastBM
fastAnc(tree,x) # estimate states
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
