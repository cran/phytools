\name{phenogram}
\alias{phenogram}
\title{Plot phenogram (traitgram)}
\usage{
phenogram(tree, x, fsize=1.0, ftype="reg", colors=NULL, axes=list())
}
\arguments{
  \item{tree}{an object of class "phylo", with or without a mapped discrete character.}
  \item{x}{a vector containing the states at the tips \emph{or} the states at all the tips and the internal nodes of the tree.}
  \item{fsize}{relative font size for tip labels.}
  \item{ftype}{font type - options are \code{"reg"}, \code{"i"} (italics), \code{"b"} (bold), or \code{"bi"} (bold-italics).}
  \item{colors}{colors for plotting the mapped character (if available) in \code{tree}.}
  \item{axes}{list of axis dimensions. Items are \code{time} and \code{trait}.}
}
\description{
  Function plots a traitgram, that is, a projection of the phylogenetic tree in a space defined by phenotype (on the \emph{y} axis) and time (on the \emph{x}). If a discrete character is mapped on the tree this will also be plotted.
}
\value{
  plots a traitgram, optionally with a mapped discrete character.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
  tree<-pbtree(n=20,scale=2)
  x<-fastBM(tree)
  phenogram(tree,x)
  # or, simulate a discrete character history
  tree<-sim.history(tree,Q=matrix(c(-1,1,1,-1),2,2),anc="1")
  # simulate in which the rate depends on the state
  x<-sim.rates(tree,c(1,10))
  phenogram(tree,x)
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
