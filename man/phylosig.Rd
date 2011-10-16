\name{phylosig}
\alias{phylosig}
\title{Compute phylogenetic signal with two methods}
\usage{
phylosig(tree,x,method="K",test=FALSE,nsim=1000)
}
\arguments{
  \item{tree}{a phylogenetic tree in \code{"phylo"} format.}
  \item{x}{vector containing values for a single continuously distributed trait.}
  \item{method}{method to compute signal: can be \code{"K"} or \code{"lambda"}.}
  \item{test}{logical indicating whether or not to conduct a hypothesis test of \code{"K"} or \code{"lambda"}.}
  \item{nsim}{for \code{method="K"}, number of simulations in randomization test.}
}
\description{
  This function computes phylogenetic signal using two different methods.  It can also conduct the hypothesis tests for significant phylogenetic signal.
}
\details{
  Lambda optimization is performed using \code{\link{optimize}} with the range of lambda set between 0 and the theoretical upper limit of lambda (determined by the relative height of the most recent internal node on the tree).
}
\value{
  if \code{(method="K")}, a list with the following components:
  \item{K}{value of the K-statistic.}
  \item{P}{optionally, the P-value from the randomization test.}
  if \code{(method="lambda")}, a list with the following components:
  \item{lambda}{fitted value of lambda.}
  \item{logL}{log-likelihood.}
  \item{P}{P-value of the likelihood ratio test.}
}
\references{
  Pagel, M. 1999. Inferring the historical patterns of biological evolution. \emph{Nature}, \bold{401}, 877--884.

  Blomberg, S. P., T. Garland Jr., A. R. Ives. 2003. Testing for phylogenetic signal in comparative data: Behavioral traits are more labile. \emph{Evolution}, \bold{57},
  717--745.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
}
\examples{
  tree<-rbdtree(1,0,Tmax=4)
  x<-fastBM(tree)
  phylosig(tree,x,method="lambda",test=TRUE)
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{simulation}
