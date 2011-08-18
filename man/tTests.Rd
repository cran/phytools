\name{tTests}
\alias{tTests}
\title{Generate all bi- and multifurcating unrooted trees}
\usage{
tTests(x,y)
}
\arguments{
  \item{x}{a vector containing the groups.}
  \item{y}{a vector containing the response variable (continuously valued).}
}
\description{
  This function computes all pairwise t-values.  Used internally by \code{\link{phylANOVA}}.
}
\details{
  Uses a little bit of code from \code{\link{pairwise.t.test}}.
}
\value{
  \item{T}{matrix of t-values.}
}
\references{
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{anova}}, \code{\link{pairwise.t.test}}, \code{\link{phylANOVA}}
}
\examples{
}
\keyword{phylogenetics}
