\name{writeNexus}
\alias{writeNexus}
\title{Write a stochastic character mapped tree to file}
\usage{
	writeNexus(tree, file="")
}
\arguments{
 	\item{tree}{object of class \code{"phylo"} or \code{"multiPhylo"}.}
	\item{file}{file name for output.}
}
\description{
	This function writes one or multiple phylogenetic trees to file in NEXUS format. Redundant with \code{ape::\link{write.nexus}}.
}
\value{
	Trees written to file.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{write.simmap}}, \code{\link{write.nexus}}
}
\keyword{phylogenetics}
\keyword{input/output}
