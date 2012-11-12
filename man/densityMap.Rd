\name{densityMap}
\alias{densityMap}
\title{Plot posterior density of stochastic mapping on a tree}
\usage{
	densityMap(trees, res=100, fsize=NULL, ftype=NULL, lwd=3, check=FALSE,
     legend=NULL, outline=FALSE)
}
\arguments{
 	\item{trees}{set of phylogenetic trees in a modified \code{"multiPhylo"} object. Values for a discrete character are mapped on the tree. See \code{\link{make.simmap}} and \code{\link{read.simmap}} for details.}
	\item{res}{resolution for gradient plotting. Larger numbers indicate a finer (smoother) gradient.}
	\item{fsize}{relative font size - can be a vector with the second element giving the font size for the legend.}
	\item{ftype}{font type - see options in \code{\link{plotSimmap}}. As with \code{fsize}, can be a vector with the second element giving font type for the legend.}
	\item{lwd}{line width for branches.}
	\item{check}{check to make sure that the topology and branch lengths of all phylogenies in \code{trees} are equal.}
	\item{legend}{if \code{FALSE} no legend is plotted; if a numeric value, it gives the length of the legend in units of branch length. Default is 0.5 times the total tree length.}
	\item{outline}{logical value indicating whether or not to outline the branches of the tree in black.}
}
\description{
	Function plots a tree with the posterior density for a mapped character from stochastic character mapping on the tree. Since the mapped value is the probability of being in state "1", only binary [0,1] characters are allowed.
}
\value{
	plots a tree.
}
\references{
	Revell, L. J. \emph{Submitted.} Two new graphical methods for mapping trait evolution on phylogenies. 

	Huelsenbeck, J. P., R. Neilsen, and J. P. Bollback. 2003. Stochastic mapping of morphological characters. \emph{Systematic Biology}, \bold{52}, 131-138. 

	Bollback, J. P. 2006. Stochastic character mapping of discrete traits on phylogenies. \emph{BMC Bioinformatics}, \bold{7}, 88.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{plotSimmap}}, \code{\link{read.simmap}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
