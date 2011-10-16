\name{write.simmap}
\alias{write.simmap}
\title{Simulate stochastic character map on a phylogenetic tree}
\usage{
write.simmap(tree,file=NULL,append=FALSE,map.order=NULL)
}
\arguments{
  \item{tree}{a phylogenetic tree as a modified object of class \code{"phylo"}. See \code{\link{make.simmap}} and \code{\link{read.simmap}}.}
  \item{file}{an optional filename.}
  \item{append}{a logical value indicating whether to append to file.}
  \item{map.order}{a optional value specifying whether to write the map in left-to-right or right-to-left order.  Acceptable values are "left-to-right" or "right-to-left" or some abbreviation of either.  If not provided, \code{write.simmap} will use \code{attr(tree,"map.order")} if available.}
}
\description{
	This function writes stochastic character mapped trees to file using the Newick style format of SIMMAP v1.0 (Bollback 2006). Note, can only write one tree at a time to file (hence the \code{append} option).
}
\details{
}
\value{
  a file or string (if \code{file=NULL}).
}
\references{
  Huelsenbeck, J. P., R. Neilsen, and J. P. Bollback. 2003. Stochastic mapping of morphological characters. \emph{Systematic Biology}, \bold{52}, 131-138. 

  Bollback, J. P. 2006. Stochastic character mapping of discrete traits on phylogenies. \emph{BMC Bioinformatics}, \bold{7}, 88.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{plotSimmap}}
}
\examples{
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{input/output}
