\name{midpoint.root}
\alias{midpoint.root}
\title{Midpoint root a phylogeny}
\usage{
midpoint.root(tree)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
}
\description{
	This function midpoint roots a rooted or unrooted tree.
}
\details{
	Midpoint rooting involves locating the midpoint of the longest path between any two tips and putting the root in that location. This function performs the same operation as \code{midpoint} in the phangorn package, but uses no phangorn code internally.
}
\value{
	A phylogenetic tree in \code{"phylo"} format.
}
\references{
	Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289--290.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{reroot}}, \code{\link{root}}
}
\keyword{phylogenetics}
\keyword{utilities}
