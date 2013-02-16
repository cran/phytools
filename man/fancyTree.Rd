\name{fancyTree}
\alias{fancyTree}
\title{Plots special types of phylogenetic trees}
\usage{
fancyTree(tree, type=c("extinction","traitgram3d","droptip","xkcd","densitymap","contmap"),
     ..., control=list())
}
\arguments{
  \item{tree}{a phylogenetic tree in \code{"phylo"} format.}
  \item{type}{the type of special plot to create. See Description.}
  \item{...}{arguments to be passed to different methods.}
  \item{control}{a list of control parameters, depending on \code{type}.}
}
\description{
  This function plots different types of phylogenetic trees. If \code{type="extinction"} (or any unambiguous abbreviation) the function will plot a tree in which branches preceding the MRCA of all extant taxa and branches leading only to extnct lineages are plotted with dashed red lines. If \code{type="traitgram3d"} the function will plot a three dimensional traitgram (that is, a projection of the tree into three dimensional morphospace where two dimensions are the phenotypic trait and the third axis is time since the root). In this case, the additional argument \code{X}, a matrix containing the tip values of all species (with species IDs as row names) should be supplied. Optionally, the user can also supply the matrix \code{A}, which contains the ancestral states in the tree with rows labeled by node number. If \code{type="droptip"} the function will create a two panel figure in which the first panel is the tree with lineages to be pruned highlighted; and the second panel is the pruned tree. In this case, the additional argument \code{tip}, the tip name or vector of tip names to be dropped, must be supplied. If \code{type="xkcd"}, an xkcd.com style tree will be plotted. (See \code{\link{xkcdTree}} for additional optional arguments if \code{type="xkcd"}.) If \code{type="densitymap"}, a posterior probability density "heat-map" is created based on a set of trees in a \code{"multiPhylo"} object containing a binary [0,1] mapped character. (See \code{\link{densityMap}} for additional optional arguments if \code{type="densitymap"}.) Finally, if \code{type="contmap"}, reconstructed continuous trait evolution is mapped on the tree. Again, see \code{\link{contMap}} for additional arguments if \code{type="contmap"}. Presently only \code{type="traitgram3d"} uses the list \code{control} which can be supplied the same set of control parameters as \code{\link{phylomorphospace3d}}, as well as the control parameter \code{maxit} which will be passed to \code{\link{anc.ML}}.
}
\value{
  This function plots different types of phylogenetic trees. For \code{type="droptip"} the function also returns the pruned tree.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{contMap}}, \code{\link{densityMap}}, \code{\link{drop.tip}}, \code{\link{phenogram}}, \code{\link{phylomorphospace3d}}, \code{\link{plot.phylo}}, \code{\link{plotSimmap}}, \code{\link{xkcdTree}}
}
\examples{
	# plot tree with extinction
	set.seed(10)
	tree<-rlineage(1,0.5,4)
	fancyTree(tree,type="extinction")

	\dontrun{
		# plot 3D traitgram
		tree<-pbtree(n=50,scale=10)
		Y<-sim.corrs(tree,vcv=matrix(c(1,0.75,0.75,1),2,2))
		fancyTree(tree,type="traitgram3d",X=Y,control=list(spin=FALSE))
	
		# plot with internal nodes from simulation
		Y<-sim.corrs(tree,vcv=matrix(c(1,0.75,0.75,1),2,2),internal=TRUE)
		B<-Y[length(tree$tip)+1:tree$Nnode,]; Y<-Y[1:length(tree$tip),]
		fancyTree(tree,type="traitgram3d",X=Y,A=B,control=list(simple.axes=TRUE,spin=FALSE))
	}
	
	# plot with dropped tips
	tree<-pbtree(n=30)
	tips<-sample(tree$tip.label)[1:10]
	pruned<-fancyTree(tree,type="droptip",tip=tips)
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{utilities}
