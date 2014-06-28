\name{locate.yeti}
\alias{locate.yeti}
\title{Locate a cryptic, recently extinct, or missing taxon on a tree}
\usage{
locate.yeti(tree, X, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{X}{a matrix with continuous character data.}
	\item{...}{optional arguments including: \code{method} (\code{"heuristic"} or \code{"exhaustive"}); \code{constraint}, a vector containing the daughter node numbers from \code{tree$edge} for each edge to try.}
}
\description{
	This function uses ML to place a recently extinct, cryptic, or missing taxon on an ultrametric tree.
}
\value{
	Optimized tree as an object of class \code{"phylo"}.
}
\references{
	Felsenstein, J. (1981) Maximum likelihood estimation of evolutionary trees from continuous characters. \emph{American Journal of Human Genetics}, 25, 471-492.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{inference}
\keyword{maximum likelihood}
