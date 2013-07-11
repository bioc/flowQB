\name{rPS}
\alias{rPS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{rPS
%%  ~~function to do ... ~~
}
\description{
This function uses Robust Statistics  to extract the statistics for the peaks found by Kmeans. 
This function process Kmeans assuming nClusters peaks to be defined and output their Robust Statistics. 
}
\usage{
rPS(TSSINGLETS, OSSINGLETS, nClusters)
}
\arguments{
  \item{TSSINGLETS}{
TSSINGLETS table having the LOGICLE events VALUES of the selected marker including SSC in the first Column. 
}
\item{OSSINGLETS}{
OSSINGLETS table having the ORIGINAL events VALUES of the selected marker including SSC in the first Column. 
}
  \item{nClusters}{
Number of the peaks to define.
}
}
\details{
Use logicle transformation to transform the data. Second identify bead singlets. Select the marker with SSC in the first Column. 
Store the 2D data in TSSINGLETS and the no-transformed 2D data in OSSINGLETS.
% 
}
\value{
The robust statistics and the number of events of each peak are returned. 
}
\references{
Faysal El Khettabi et al. 2013, Automated Quadratic Characterization of Flow Cytometer Instrument Sensitivity, to be submitted. 
}
\author{
Faysal el Khettabi 
}
\note{
This function is a 2D processing of selected marker with SSC. A general function is given in the function MultilevelrPS.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See MultilevelrPS function
}
\examples{
# see the  vignette flowQB_Sweave_NIH2_.Rnw for a detailed example.
}

