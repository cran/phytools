\name{fastBM}
\alias{fastBM}
\title{Fast Brownian simulation}
\usage{
fastBM(tree,a=0,mu=0,sig2=1,bounds=c(-Inf,Inf),internal=FALSE,nsim=1)
}
\arguments{
  \item{tree}{is a phylogenetic tree in \code{"phylo"} format.}
  \item{a}{a value for ancestral state at the root node.}
  \item{mu}{an optional value for the mean of random normal changes along branches of the tree - can be used to simulate a trend if \code{mu!=0}.}
  \item{sig2}{instantaneous variance of the BM process.}
  \item{bounds}{a vector with the lower and upper bounds (respectively) for bounded Brownian simulation - by default simulation is unbounded.}
  \item{internal}{logical value indicating whether or not to return states for internal nodes.}
  \item{nsim}{number of simulations.}
}
\description{
  This function conducts fast BM simulation on a tree.
}
\details{
}
\value{
  a vector (for \code{nsim=1}) or matrix containing the tip states for the \code{n} species in the tree, and (optionally) the ancestral states for internal nodes.
}
\references{
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{branching.diffusion}}
}
\examples{
}
\keyword{phylogenetics}
\keyword{simulation}
