\name{phylANOVA}
\alias{phylANOVA}
\title{Phylogenetic ANOVA and post-hoc tests}
\usage{
phylANOVA(tree, x, y, nsim=1000, posthoc=TRUE, p.adj="holm")
}
\arguments{
	\item{tree}{a phylogenetic tree in \code{"phylo"} format.}
	\item{x}{a vector containing the groups.}
	\item{y}{a vector containing the response variable (continuously valued).}
	\item{nsim}{an integer specifying the number of simulations (including the observed data).}
	\item{posthoc}{a logical value indicating whether or not to conduct posthoc tests to compare the mean among groups.}
	\item{p.adj}{method to adjust P-values for the posthoc tests to account for multiple testing. Options same as \code{\link{p.adjust}}.}
}
\description{
	Simulation based phylogenetic ANOVA following Garland et al. (1993), with post-hoc tests.
}
\details{
	This function performs the simulation-based phylogenetic ANOVA of Garland et al. (1993) and (optionally) conducts all posthoc comparisons of means among groups (also obtaining the P-values by phylogenetic simulation).

	This function uses a little bit of code from both \code{phy.anova} in the \pkg{geiger} package and \code{\link{pairwise.t.test}}.
}
\value{
	An object of class \code{"phylANOVA"} containing the following elements:
	\item{F}{F from observed data.}
	\item{Pf}{P-value for F from simulation.}
	\item{T}{matrix of t-values.}
	\item{Pt}{matrix of multiple test corrected P-values from posthoc t-tests.}
}
\references{
	Garland, T., Jr., A. W. Dickerman, C. M. Janis, & J. A. Jones (1993) Phylogenetic analysis of covariance by computer simulation. \emph{Systematic Biology}, \bold{42}, 265-292.

	Harmon, L. J., J. T. Weir, C. D. Brock, R. E. Glor, W. Challenger (2008) GEIGER: investigating evolutionary radiations. \emph{Bioinformatics}, 24, 129-131.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{anova}}, \code{\link{pairwise.t.test}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{statistics}
\keyword{least squares}
\keyword{simulation}
