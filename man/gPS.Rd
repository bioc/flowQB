\name{gPS}
\alias{gPS}
\title{gPS
Process Kmeans and fit the data to a Gaussian distribution using MASS package.
}
\description{
This function uses Gaussian distribution fitting  to extract the statistics for the peaks found by Kmeans. 
This function is usefull when a marker is associated to SSC to form a 2D processing using Kmeans.
This function process Kmeans assuming nClusters peaks to be detected and output their statistics with the associated number of events. 
}
\usage{
gPS(TSSINGLETS, OSSINGLETS, nClusters)
}
\arguments{
  \item{TSSINGLETS}{
TSSINGLETS table having the LOGICLE events VALUES of the selected marker including SSC in the first Column. 
}
\item{OSSINGLETS}{
OSSINGLETS table having the ORIGINAL events VALUES of the selected marker including SSC in the first Column. 
}
  \item{nClusters}{
Number of the peaks to detect.
}
}
\details{
Use logicle transformation to transform the data. Second identify bead singlets. Select the marker with SSC in the first Column. 
Store the 2D data in TSSINGLETS and the no-transformed 2D data in OSSINGLETS.
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
This function is a 2D processing of selected marker with SSC. A general function is given in the function MultilevelgPS.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See MultilevelgPS function
}
\examples{
# The vignettes $AdvancedflowQBNIH2.Rnw$ and $AdvancedflowQBNIH3.Rnw$ have a detailed illustration how to use this function.
}
