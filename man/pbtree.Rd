\name{pbtree}
\alias{pbtree}
\title{Simulate pure-birth stochastic tree or trees}
\usage{
pbtree(b=1,n=NULL,t=NULL,scale=NULL,nsim=1,ape=TRUE)
}
\arguments{
  \item{b}{birth rate or speciation rate.}
  \item{n}{desired number of species.}
  \item{t}{time stopping criterion (not yet implemented).}
  \item{scale}{if set, rescales tree to have total length \code{scale}.}
  \item{nsim}{number of simulated trees to return.}
  \item{ape}{return nodes in a 'ape' compatible ordering.}
}
\description{
  This function simulates stochastic pure-birth (Yule) trees. Note that if \code{ape=FALSE}, then the function will run faster, and the tree is theoretically compatible with the ape \code{'phylo'} standard; however some downstream errors with functions such as \code{link{bind.tree}} have been observed.
}
\value{
  a tree or list of trees.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
   # simulate a pure-birth tree with 1000 tips, scaled to a length of 1.0
   tree<-pbtree(n=1000,scale=1)
}
\keyword{phylogenetics}
\keyword{simulation}
