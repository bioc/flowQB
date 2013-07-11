\name{lrMESF}
\alias{lrMESF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ lrMESF
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function lrMESF uses the Means and Variance in terms of the MESF values to conduct a linear regression. The obtained coefficients are used to calculate the Q and B values.
}
\usage{
lrMESF(MESFmeanssd2, Peak1, Peak2)
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
The detector efficiency (Q) and optical background (B) are derived from the obtained linear regression coefficients when the SDs are corrected for illumination.
}
\value{
Output linear regression detector efficiency (Q) and optical background (B).
}
\author{Faysal El Khettabi \\
Terry Fox Laboratory \\
British Columbia Cancer Agency \\
Vancouver, BC, Canada \\
E-mail: fkhettabi@bccrc.ca 
}

\references{ 
Faysal El Khettabi et al. 2013, Automated Quadratic Characterization of Flow Cytometer Instrument Sensitivity, to be submitted. }
\keyword{ package}
\seealso{ 
J. Wood, Fundamental Flow Cytometer Properties Governing Sensitivity and Resolution, 
Cytometry 33, (1998), p.~ 260 - 6.

E. Chase and R. Hoffman, Resolution of Dimly Fluorescent Particles: a Practical Measure of Fluorescence Sensitivity,
Cytometry 33 (1998), p.~ 267-279.

R. Hoffman and J. Wood, Characterization of Flow Cytometer Instrument Sensitivity, Current Protocols in Cytometry, Chapter 1: Unit 1.20 (2007).

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

# Usage of the function ReadDD	
t2D=ReadDD(File,1,2,P,5,3)
summary(t2D)

# Usage of the function KmeansMeanSD	
MFIMeansSDs=KmeansMeanSD(t2D,8,500,200,1)
# Function KmeansMeanSD returns the MFI Means and SDs of the 8 clusters.
MFIMeansSDs

# For MESF calculation, the constant conversion between MFI and MESF is set to:
p=357217.00/7102

# MFIs are converted to MESFs and SDs are not corrected as we set,
IllCorrCV=MFIMeansSDs[8,2]/MFIMeansSDs[8,1] #SDs are corrected using the beads in cluster 8.

# Usage of the function MFI2MESF	
MFI2MSEF=MFI2MESF(MFIMeansSDs,p,IllCorrCV)
MFI2MSEF
# Note that MESF: MESF Mean and MESFV: MESF Variance (SD^2).

# Linear regression 
# Peaks associated to cluster  3 to  cluster 5 are used 
# to compute the linear regression coefficients.

# Usage of the function lrMESF
LQB=lrMESF(MFI2MSEF,3,5)
print("Linear QB")
OV=c(Q=as.numeric(LQB[1]), B=as.numeric(LQB[2]), Rsquared=as.numeric(round(LQB[3], 2)))
print(OV)
# Note  c1 =1/Q, c0= B/Q 



}
