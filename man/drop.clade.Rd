\name{drop.clade}
\alias{drop.clade}
\title{Drop clade from a tree}
\usage{
drop.clade(tree,tip)
}
\arguments{
  \item{tree}{phylogeny.}
  \item{tip}{set of tips in a clade.}
}
\description{
  Mostly internal function for \code{\link{posterior.evolrate}}; function drops the clade containing the species in \code{tip}.
}
\details{
  Probably should not use unless you know what you're doing.
}
\value{
  a tree.
}
\references{
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
}
\examples{
}
\keyword{phylogenetics}
\keyword{utilities}
