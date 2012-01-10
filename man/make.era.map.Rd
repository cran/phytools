\name{make.era.map}
\alias{make.era.map}
\title{Create "era" map on a phylogenetic tree}
\usage{
make.era.map(tree,limits)
}
\arguments{
  \item{tree}{a phylogenetic tree as an object of class \code{"phylo"}.}
  \item{limits}{a vector containing the temporal limits, in time since the root node of the tree, for the mappings.}
}
\description{
	This function creates a temporal map on the tree based on \code{limits} provided by the user.
}
\details{
}
\value{
  a modified phylogenetic tree of class \code{"phylo"} with the following additional elements:
  \item{maps}{a list of named vectors containing the times spent in each state on each branch, in the order in which they occur.}
  \item{mapped.edge}{a matrix containing the total time spent in each state along each edge of the tree.}
}
\references{
  Paradis, E., J. Claude, and K. Strimmer. 2004. APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289--290.

}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{plotSimmap}}
}
\examples{
}
\keyword{phylogenetics}
\keyword{comparative method}

