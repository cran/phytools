\name{bmPlot}
\alias{bmPlot}
\title{Simulates and visualizes discrete-time Brownian evolution on a phylogeny}
\usage{
bmPlot(tree, type="BM", anc=0, sig2=1/1000, ngen=1000, ...)
}
\arguments{
  \item{tree}{a phylogenetic tree in \code{"phylo"} format.}
  \item{type}{the type of plot to create. See Description.}
  \item{anc}{the ancestral value for the root node.}
  \item{sig2}{the BM rate (variance of the Brownian evolution process).}
  \item{ngen}{number of generations for the simulation: will rescale the tree to this total length.}
  \item{...}{arguments to be passed to different methods.}
}
\description{
  This function conducts discrete-time Brownian motion simulation on an input tree, plots the outcome, and returns the tip and internal node states to the user as a named vector. The function will first rescale and round the branch lengths to integer length. If integer branch lengths are provided, the user should also set \code{ngen=max(nodeHeights(tree))}. For \code{type="threshold"} the visualization is of the threshold model (Felsenstein 2012), in which the evolving character is liability and the segements of evolution are colored by their value for the threshold trait. If \code{type="threshold"} is used, the function requires at least one addition input: \code{thresholds}, a vector containing the ordered thresholds between states. The user can also provide the colors for plotting in \code{cols}. Note that one more color than threshold should be provided as one threshold implies two states; two thresholds, three states; etc. If no value for \code{cols} is provided, the function will recycle a set of four colors up to the number of times required by \code{thresholds}.
}
\value{
  This function conducts and plots discrete time Brownian simulation and returns a vector containing the simulated states at internal nodes and tips of the tree.
}
\references{
  Felsenstein, J. 2012. A comparative method for both discrete and continuous characters using the threshold model. \emph{American Naturalist}, \bold{179}, 145-156.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{fastBM}}, \code{\link{phenogram}}, \code{\link{threshBayes}}
}
\examples{
	# plot BM simulation on 20 taxon tree
	tree<-pbtree(n=20)
	x<-bmPlot(tree)
	
	# plot simulation of a threshold character
	tree<-pbtree(n=20)
	x<-bmPlot(tree,type="threshold",thresholds=c(0,1,2))

}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{simulation}
