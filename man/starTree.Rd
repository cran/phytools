\name{starTree}
\alias{starTree}
\title{Create star phylogeny}
\usage{
starTree(species,branch.lengths=NULL)
}
\arguments{
  \item{species}{a list of species.}
  \item{branch.lengths}{an optional list of branch lengths in the same order as \code{species}.}
}
\description{
  This function creates a star phylogeny.
}
\details{
  Creates a star phylogeny with (optionally) user specified branch lengths.
}
\value{
  a phylogenetic tree in \code{"phylo"} format.
}
\references{
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{stree}}
}
\examples{
}
\keyword{phylogenetics}
