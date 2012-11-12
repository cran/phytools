\name{branching.diffusion}
\alias{branching.diffusion}
\title{Animation of branching random diffusion}
\usage{
branching.diffusion(sig2=1, b=0.0023, time.stop=1000, ylim=NULL, 
     smooth=TRUE, pause=0.02, record=NULL, path=NULL)
}
\arguments{
  \item{sig2}{variance of BM process.}
  \item{b}{birthrate for branching process.}
  \item{time.stop}{number of generations to run.}
  \item{ylim}{y limits (for plotting).}
  \item{smooth}{if TRUE, smoother plotting (but slower).}
  \item{pause}{pause (in s) between generations.}
  \item{record}{filename for video file output (no video if \code{NULL}).}
  \item{path}{path to file for video rendering (by default is \code{C:/Program Files/ffmpeg/bin/ffmpeg.exe}).}
}
\description{
  This function creates an animation of branching random diffusion (i.e., BM with speciation).
}
\value{
  None.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{fastBM}}
}
\keyword{phylogenetics}
\keyword{animation}
\keyword{plotting}
