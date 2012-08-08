\name{pgls.Ives}
\alias{pgls.Ives}
\title{Phylogenetic regression with intraspecific sampling error}
\usage{
pgls.Ives(tree,X,y,Vx,Vy,Cxy)
}
\arguments{
	\item{tree}{a phylogenetic tree in \code{"phylo"} format.}
	\item{X}{a named vector containing a \emph{single} independent variable (multiple independent variables to be added).}
	\item{y}{vector the dependent variable.}
	\item{Vx}{sampling variances for \code{X}.}
	\item{Vy}{sampling variances for \code{y}.}
	\item{Cxy}{sampling covariances between \code{X} and \code{y}.}
}
\description{
	This function fits the phylogenetic regression model with within-species sampling error following Ives et al. (2007).
}
\details{
	Presently only the bivariate regression model is implemented. Uses \code{optim(...,method="L-BFGS-B")}.
}
\value{
	a list with the following components:
	\item{beta}{a vector or matrix of regression coefficients.}
	\item{sig2x}{fitted BM rate for \code{X}.}
	\item{sig2y}{fitted BM rate for \code{y}.}
	\item{a}{fitted ancestral states for \code{X} and \code{y}.}
	\item{logL}{log-likelihood.}
}
\references{
	Ives, A. R., P. E. Midford, and T. Garland Jr. 2007. Within-species measurement error in phylogenetic comparative methods. \emph{Systematic Biology}, \bold{56},
  252--270.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownie.lite}}, \code{\link{phylosig}}, \code{\link{phyl.resid}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{statistics}
\keyword{least squares}
\keyword{maximum likelihood}
