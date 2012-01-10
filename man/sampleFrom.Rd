\name{sampleFrom}
\alias{sampleFrom}
\title{Sample from a set of distributions}
\usage{
sampleFrom(xbar=0,xvar=1,n=1,randn=NULL,type="norm")
}
\arguments{
  \item{xbar}{a named vector of means.}
  \item{xvar}{a named vector of variances.}
  \item{n}{a vector containing the sample sizes of each species.}
  \item{randn}{a range of sample sizes are to be random.}
  \item{type}{\code{"norm"} is the only distribution implemented so far.}
}
\description{
  Function samples from a set of normal distributions with parameters given in \code{xbar} and \code{xvar}.
}
\details{
}
\value{
  a vector, with labels.
}
\references{

}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
}
\examples{

}
\keyword{phylogenetics}
\keyword{statistics}