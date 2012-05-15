\name{DiscriminantExamination}
\alias{DiscriminantExamination}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ DiscriminantExamination
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
DiscriminantExamination uses the values in the output of the function qrMESF to estimate 
the discriminant of the resulting quadratic equation and can be used as an additional interpretation tool to aid in understanding cytometer
sensitivity.

}
\usage{
DiscriminantExamination(Q, B, sigmaS2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{
%%     ~~Describe \code{Q} here~~
Q: the first value in the output of the function qrMESF.
}
  \item{B}{
%%     ~~Describe \code{B} here~~
B: the second  value in the output of the function qrMESF.
}
\item{sigmaS2}{
%%     ~~Describe \code{Peak2} here~~
sigmaS2: the fourth  value in the output of the function qrMESF.
}
}


\details{
It examines the discriminant, Delta = c1^2 -4*c0*c2, of the quadratic equation, c0+c1*MESF+c2*MESF^2.

There are two possible scenarios: Delta >= 0 or Delta  < 0. The coefficients (c0,c1,c2) are derived from the output of the function qrMESF.
}
\value{
Output c0, c1, c2 and Delta values. 

If Delta >=  0, the larger the variation product of (sigmaE2) and (sigmaS2), the lower the upper bound on the detection efficiency Q.

If Delta < 0, the lower the variation product of (sigmaE2) and (sigmaS2), the greater the upper bound on the detection efficiency Q. 
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
The function qrMESF in flowQB.
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

# Discriminant of the Quadratic Equation
Coefs <- DiscriminantExamination(as.real(QQB[1]),as.real(QQB[2]),as.real(QQB[4]))

Delta=Coefs[4]

if(Delta >= 0) 
{
cat(paste("The sign of the discriminant is positive with the value",round(Delta,2),"\n"))
cat("The larger the variation product of (sigmaE2) and (sigmaS2), \n")
cat("the lower the upper bound on the detection efficiency Q.\n")
} 

if(Delta < 0) 
{
cat(paste("The sign of the discriminant is negative with the value",round(Delta,2),"\n"))
cat("The lower the variation product of (sigmaE2) and (sigmaS2), \n")
cat("the greater the upper bound on the detection efficiency Q. \n")
}



}
