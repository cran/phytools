\name{reroot}
\alias{reroot}
\title{Re-root a tree along an edge}
\usage{
reroot(tree,node.number,position)
}
\arguments{
  \item{tree}{a phylogenetic tree in \code{"phylo"} format.}
  \item{node.number}{ number of the node descending from the target branch in \code{tree$edge} - this can also be a tip in which case the node number is the index number of the tip in \code{tree$tip.label}.}
  \item{position}{position along the target edge at which to re-root the tree.}.
}
\description{
  This function re-roots a phylogenetic tree at an arbitrary position along an edge.
}
\details{
  This function essentially acts as a mildly sophisticated wrapper for \code{\link{root}} in the "ape" package.
}
\value{
  a phylogenetic tree in \code{"phylo"} format.
}
\references{
  Paradis, E., J. Claude, and K. Strimmer. 2004. APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289--290.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{root}}
}
\keyword{phylogenetics}
\keyword{utilities}
