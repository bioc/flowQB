\name{MFI2MESF}
\alias{MFI2MESF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MFI2MESF
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
MFI2MESF converts the MFI means to MESF means.
}
\usage{
MFI2MESF(MeansSDs, p, IllCorrCV)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MeansSDs}{
2D array, first column MFI means , second column the associated SDs.
}
  \item{p}{
Constant Conversion between MFI and MESF
}
  \item{IllCorrCV}{
Value to be used to correct the SDs.
}
}
\details{
The output will be used as input in the functions lrMESF or qrMESF to calculate the detector efficiency (Q) and optical background (B).
}
\value{
Convert MFI to MSEF, the output is a 2D array having the MESF Means (MESF) and the associated MESF Variances (MESFV).
}
\references{
Faysal El Khettabi et al. 2012, Flow Cytometer Sensitivity: A Quadratic Model, to be submitted. 

A. Gaigalas and L. Wang, Approaches to Quantitation in Flow Cytometry, in Standardization and Quality Assurance in Fluorescence Measurements II Springer Series on Fluorescence (2008), Volume 6, Part D, 371-398.
}
\author{Faysal El Khettabi \\
Terry Fox Laboratory \\
British Columbia Cancer Agency \\
Vancouver, BC, Canada \\
E-mail: fkhettabi@bccrc.ca 
}

\seealso{ 
The functions lrMESF and  qrMESF in flowQB.
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




}
