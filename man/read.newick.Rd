\name{read.newick}
\alias{read.newick}
\title{Read Newick style tree}
\usage{
read.newick(file="",text)
}
\arguments{
  \item{file}{name of text file with single Newick style tree.}
  \item{text}{character string containing tree.}
}
\description{
  This function reads a Newick style tree from file.
}
\details{
  Note, this function does not do anything that \code{\link{read.tree}} doesn't do better - and is included mostly for development purposes.
}
\value{
  a phylogenetic tree of class \code{"phylo"}
}
\references{
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{read.tree}}, \code{\link{read.nexus}}
}
\examples{
tree<-"((Human,Chimp),Gorilla),Monkey);"
phy<-read.newick(text=tree)
}
\keyword{phylogenetics}
\keyword{input/output}
