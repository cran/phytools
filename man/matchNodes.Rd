\name{matchNodes}
\alias{matchNodes}
\title{Matches nodes between two trees}
\usage{
	matchNodes(tr1, tr2, method=c("descendants","distances"))
}
\arguments{
	\item{tr1}{first tree.}
	\item{tr2}{second tree.}
	\item{method}{method to use to match nodes between trees. \code{"descendants"} uses the tip species descended from each node; \code{"distances"} uses the distances from the nodes to the tips.}
}
\description{
	This function returns a matrix in which the first column contains all the internal nodes for \code{tr1} and the second column contains the matching nodes from \code{tr2}, inasmuch as they can be identified. A pair of matching nodes are defined by sharing all descendant leaves in common.
}
\details{
	Primarily to be used internally by \code{\link{fastAnc}}.
}
\value{
	a matrix.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
