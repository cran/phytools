\name{mrp.supertree}
\alias{mrp.supertree}
\title{Matrix representation parsimony supertree estimation}
\usage{
mrp.supertree(phy, weights=NULL)
}
\arguments{
  \item{phy}{an object of class \code{"multiPhylo"} (i.e., a list of trees).}
  \item{weights}{not used.}
}
\description{
  This function estimates the MRP supertree from a set of trees.
}
\details{
  Function uses \code{pratchet} from the "phangorn" package (Schliep 2011) and \code{prop.part} from the "ape" package (Paradis et al. 2004).
}
\value{
  a \code{"phylo"} or \code{"multiPhylo"} object that is the MP or set of MP MRP trees.
}
\references{
  Felsenstein, J. 2004. \emph{Inferring Phylogenies}. Sinauer.

  Paradis, E., J. Claude, and K. Strimmer. 2004. APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289--290.

  Schliep, K. P. 2011. phangorn: phylogenetic analysis in R. \emph{Bioinformatics}, \bold{27}, 592-593.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{exhaustiveMP}}, \code{\link{optim.parsimony}}, \code{\link{pratchet}}
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{parsimony}
\keyword{supertree}
