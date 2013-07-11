\name{KmeansMedianrSD}
\alias{KmeansMedianrSD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{KmeansMedianrSD}
\description{
Function KmeansMedianrSD takes a given 2D array and generates a number of clusters and returns their Medians and robust standard deviations (rSDs).
}
\usage{
KmeansMedianrSD(transformed2Darray, nClusters, nstart, itermax, Vis)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{transformed2Darray}{
%%     ~~Describe \code{transformed2Darray} here~~
2D array.
}
  \item{nClusters}{
%%     ~~Describe \code{nClusters} here~~
Number of clusters
}
  \item{nstart}{
%%     ~~Describe \code{nstart} here~~
See R kmeans. 
}
  \item{itermax}{
%%     ~~Describe \code{itermax} here~~
See R kmeans. 
}
  \item{Vis}{
%%     ~~Describe \code{Vis} here~~
If 0 no visualization, if 1 with visualization.
}
}
\details{
See the manuscript in the references.
}
\value{
Output for each cluster the robust statistics: Medians and robust standard deviations (rSDs).
}
\references{
%% ~put references to the literature/web site here ~
Faysal El Khettabi et al. 2013, Automated Quadratic Characterization of Flow Cytometer Instrument Sensitivity, to be submitted. 
}
\author{Faysal El Khettabi \\
Terry Fox Laboratory \\
British Columbia Cancer Agency \\
Vancouver, BC, Canada \\
E-mail: fkhettabi@bccrc.ca \\
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
kmeans in R software.
}
\examples{
# see Examples.
}
