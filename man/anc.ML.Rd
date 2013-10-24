\name{anc.ML}
\alias{anc.ML}
\title{Ancestral character estimation using likelihood}
\usage{
anc.ML(tree, x, maxit=2000)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a vector of tip values for species; \code{names(x)} should be the species names.}
	\item{maxit}{an optional integer value indicating the maximum number of iterations for optimization.}
}
\description{
	This function estimates the evolutionary parameters and ancestral states for Brownian evolution using likelihood.
}
\details{
	Because this function relies on a high dimensional numerical optimization of the likelihood function, \code{\link{fastAnc}} should probably be preferred for most purposes. If using \code{\link{anc.ML}}, users should be cautious to ensure convergence. This has been ameliorated in phytools>=0.2-48 by seeding the ML optimization with the result from \code{\link{fastAnc}}.
}
\value{
	A list with the following components:
	\item{sig2}{the variance of the BM process.}
	\item{ace}{a vector with the ancestral states.}
	\item{logLik}{the log-likelihood.}
	\item{convergence}{the value of \code{$convergence} returned by \code{optim()} (0 is good).}
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ace}}, \code{\link{anc.Bayes}}, \code{\link{fastAnc}}, \code{\link{optim}}
}
\examples{
tree<-pbtree(n=50)
x<-fastBM(tree) # simulate using fastBM
anc.ML(tree,x) # fit model & estimate ancestral states
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
