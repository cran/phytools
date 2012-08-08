\name{fancyTree}
\alias{fancyTree}
\title{Plots special types of phylogenetic trees}
\usage{
fancyTree(tree,type=c("extinction","traitgram3d"),...,control=list())
}
\arguments{
  \item{tree}{a phylogenetic tree in \code{"phylo"} format.}
  \item{type}{the type of special plot to create. See Description.}
  \item{...}{arguments to be passed to different methods.}
  \item{control}{a list of control parameters, depending on \code{type}.}
}
\description{
  This function plots different types of phylogenetic trees. If \code{type="extinction"} the function will plot a tree in which branches preceding the MRCA of all extant taxa and branches leading only to extnct lineages are plotted with dashed red lines. If \code{type="traitgram3d"} the function will plot a three dimensional traitgram (that is, a projection of the tree into three dimensional morphospace where two dimensions are the phenotypic trait and the third axis is time since the root). In the latter case, the additional argument \code{X}, a matrix containing the tip values of all species (with species IDs as row names) should be supplied. Optionally, the user can also supply the matrix \code{A}, which contains the ancestral states in the tree with rows labeled by node number. Presently only \code{type="traitgram3d"} uses the list \code{control} which can be supplied the same set of control parameters as \code{\link{phylomorphospace3d}}, as well as the control parameter \code{maxit} which will be passed to \code{\link{anc.ML}}.
}
\value{
  This function plots different types of phylogenetic trees.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{phenogram}}, \code{\link{phylomorphospace3d}}, \code{\link{plot.phylo}}, \code{\link{plotSimmap}}
}
\examples{
	# plot tree with extinction
	set.seed(10)
	tree<-rlineage(1,0.5,4)
	fancyTree(tree,type="extinction")

	# plot 3D traitgram
	tree<-pbtree(n=50,scale=10)
	Y<-sim.corrs(tree,vcv=matrix(c(1,0.75,0.75,1),2,2))
	fancyTree(tree,type="traitgram3d",X=Y,control=list(spin=FALSE))
	
	# plot with internal nodes from simulation
	Y<-sim.corrs(tree,vcv=matrix(c(1,0.75,0.75,1),2,2),internal=TRUE)
	B<-Y[length(tree$tip)+1:tree$Nnode,]; Y<-Y[1:length(tree$tip),]
	fancyTree(tree,type="traitgram3d",X=Y,A=B,control=list(simple.axes=TRUE,spin=FALSE))
}
\keyword{phylogenetics}
\keyword{plotting}
