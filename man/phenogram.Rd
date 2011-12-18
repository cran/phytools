\name{phenogram}
\alias{phenogram}
\title{Plot phenogram (traitgram)}
\usage{
phenogram(tree,x,fsize=1.0,ftype="reg")
}
\arguments{
  \item{tree}{an object of class "phylo".}
  \item{x}{a vector containing the states at the tips \emph{or} the states at all the tips and the internal nodes of the tree.}
  \item{fsize}{relative font size for tip labels.}
  \item{ftype}{font type - options are \code{"reg"}, \code{"i"} (italics), \code{"b"} (bold), or \code{"bi"} (bold-italics).}
}
\description{
  Function plots a traitgram, that is, a projection of the phylogenetic tree in a space defined by phenotype (on the \emph{y} axis) and time (on the \emph{x}).
}
\details{
}
\value{
  plots a traitgram.
}
\references{

}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
}
\examples{
  tree<-rcoal(n=20)
  x<-fastBM(tree)
  phenogram(tree,x)
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
