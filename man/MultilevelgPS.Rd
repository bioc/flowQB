\name{MultilevelgPS}
\alias{MultilevelgPS}
\title{MultilevelgPS
Process multi-level Kmeans and fit the data to a Gaussian distribution using MASS package.
}
\description{
This function uses Gaussian distribution fitting  to extract the statistics for the peaks found by multi-level Kmeans.
This function is usefull when a number of markers are processed as a multidimensional space to detect the assocaited peaks for each marker in the data.  
This function process Kmeans assuming nClusters peaks to be defined and output their statistics with the associated number of events. 
}
\usage{
MultilevelgPS(TSSINGLETS, OSSINGLETS, nClusters)
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
The probability distribution fitting of the peak data events to a Gaussian distribution is used for peak statistics extraction, the fitting parameters means and standard deviations are the estimations of MFIs and SDs values with the associated number of events. 
}
\references{
Faysal El Khettabi et al. 2013, Automated Quadratic Characterization of Flow Cytometer Instrument Sensitivity, to be submitted. 
}
\author{
Faysal el Khettabi 
}
\note{
This function is a multi-level processing as the function gPS without using the SSC.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See gPS function
}
\examples{
#  The vignettes $AdvancedflowQBNIH2.Rnw$ and $AdvancedflowQBNIH3.Rnw$ have a detailed illustration how to use this function.
}

