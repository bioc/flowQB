\name{qrMESF}
\alias{qrMESF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ qrMESF
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
qrMESF uses the Means and variances in terms of the MESF values to conduct a quadratic regression
}
\usage{
qrMESF(MESFmeanssd2, Peak1, Peak2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MESFmeanssd2}{
%%     ~~Describe \code{MESFmeanssd2} here~~
MESF 2D data having means and Variances.
}
  \item{Peak1}{
%%     ~~Describe \code{Peak1} here~~
First peak to consider.
}
  \item{Peak2}{
%%     ~~Describe \code{Peak2} here~~
Last peak to consider.
}
}
\details{
The detector efficiency (Q) and optical background (B) are derived from the obtained quadratic regression coefficients.
}
\value{
Output quadratic regression detector efficiency (Q) and optical background (B).
}
\author{Faysal El Khettabi \\
Terry Fox Laboratory \\
British Columbia Cancer Agency \\
Vancouver, BC, Canada \\
E-mail: fkhettabi@bccrc.ca
}

\references{ 
Faysal El Khettabi et al. 2012, Flow Cytometer Sensitivity: A Quadratic Model, to be submitted. 
}
\keyword{ package}
\seealso{ 
lrMESF function.
}
\examples{

rm(list=ls(all=TRUE))
rm(list=ls(all=TRUE))
library("flowQB")
File= system.file("extdata","NIH.fcs",package="flowQB")

# Doublet Discriminations 
P=96

# Reading the FCS file, the Forward Scattering Area 1 and 2 
# are used to obtain singlet events with P=96%.
# The processing returns a 2D singlet events for the channel 
# of interest, index 5, with the companion channel Side Scattering, index 3.
  
library("flowQB")

# Usage of the function ReadDD	
t2D=ReadDD(File,1,2,P,5,3)
summary(t2D)

# Usage of the function KmeansMeanSD	
MFIMeansSDs=KmeansMeanSD(t2D,8,500,200,1)
# Function KmeansMeanSD returns the MFI Means and SDs of the 8 clusters.
MFIMeansSDs

# For MESF calculation, the constant conversion between MFI and MESF is set to:
p=357217.00/7102

# MFIs are converted to MESFs and SDs are not corrected as we set IllCorrCV=0.
# Usage of the function MFI2MESF	
MFI2MSEF=MFI2MESF(MFIMeansSDs,p,0)
MFI2MSEF
# Note that MESF: MESF Mean and MESFV: MESF Variance (SD^2).


# Quadratic Q and B Calculation.
# Peaks associated to cluster  3 to  cluster 6 are used to compute 
# the quadratic regression coefficients.

# Usage of the function qrMESF
QQB=qrMESF(MFI2MSEF,3,6)
print("Quadratic Q & B")
OV=c(Q=as.real(QQB[1]), B=as.real(QQB[2]), Rsquared=as.real(round(QQB[3], 2)), sigmaS2=as.real((QQB[4])))
print(OV)
# Note  c1 =1/Q, c0= B/Q and c2 = sigmaS2



}
