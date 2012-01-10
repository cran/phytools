\name{pbtree}
\alias{pbtree}
\title{Simulate pure-birth stochastic tree or trees}
\usage{
pbtree(b=1,n=NULL,t=NULL,scale=NULL,nsim=1)
}
\arguments{
  \item{b}{birth rate or speciation rate.}
  \item{n}{desired number of species.}
  \item{t}{time stopping criterion (not yet implemented).}
  \item{scale}{if set, rescales tree to have total length \code{scale}.}
  \item{nsim}{number of simulated trees to return.}
}
\description{
  This function simulates stochastic pure-birth (Yule) trees.
}
\details{
}
\value{
  a tree or list of trees.
}
\references{
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
}
\examples{
}
\keyword{phylogenetics}
\keyword{simulation}
