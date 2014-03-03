\name{vcvPhylo}
\alias{vcvPhylo}
\title{Calculates cophenetic (i.e., phylogenetic VCV) matrix}
\usage{
vcvPhylo(tree, anc.nodes=TRUE, ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{anc.nodes}{logical value indicating whether or not to include ancestral nodes.}
	\item{...}{optional arguments including \code{internal} (synonym of \code{anc.nodes}) and \code{model} (can be \code{"BM"}, \code{"OU"}, or \code{"lambda"}.}
}
\description{
	This function returns a so-called \emph{phylogenetic variance covariance matrix} (e.g., see \code{\link{vcv.phylo}}), but (optionally) including ancestral nodes and under different evolutionary models.
}
\value{
	A matrix.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{statistics}
\keyword{utilities}
