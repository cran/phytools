\name{drop.clade}
\alias{drop.clade}
\title{Drop a clade from a tree}
\usage{
drop.clade(tree, tip)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{tip}{set of tips in a clade.}
}
\description{
	Mostly internal function for \code{\link{posterior.evolrate}}; function drops the clade containing the species in \code{tip}.
}
\details{
	Probably should not use unless you know what you're doing.
}
\value{
	An object of class \code{"phylo"}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
