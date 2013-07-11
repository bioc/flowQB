\name{MultilevelrPS}
\alias{MultilevelrPS}
\title{MultilevelrPS
Process multi-level Kmeans and extract the Robust statistics.
}
\description{
This function uses Robust statistics for the peaks found by multi-level Kmeans.
This function is usefull when a number of markers are processed as a multidimensional space to detect the assocaited peaks for each marker in the data.  
This function process a multi-level Kmeans assuming nClusters peaks to be detected and output their statistics with the associated number of events. 

}
\usage{
MultilevelrPS(TSSINGLETS, OSSINGLETS, nClusters)
}
\arguments{
  \item{TSSINGLETS}{
TSSINGLETS table having the LOGICLE events VALUES of the selected markers. 
}
\item{OSSINGLETS}{
OSSINGLETS table having the ORIGINAL events VALUES of the selected markers. 
}
  \item{nClusters}{
Number of the peaks to detect.
}
}
\details{
Use logicle transformation to transform the data. Second identify bead singlets. Select the markers to be processed. 
Store the data of the selected markers in TSSINGLETS and the no-transformed data in OSSINGLETS.
% 
}
\value{
Robust statistics with the associated number of events. 
}
\references{
Faysal El Khettabi et al. 2013, Automated Quadratic Characterization of Flow Cytometer Instrument Sensitivity, to be submitted. 
}
\author{
Faysal el Khettabi 
}
\note{
This function is a multi-level processing as the function rPS without using the SSC.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See rPS function
}
\examples{
#  The vignettes $AdvancedflowQBNIH2.Rnw$ and $AdvancedflowQBNIH3.Rnw$ have a detailed illustration how to use this function.
}
