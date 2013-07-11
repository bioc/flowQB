\name{flowQB-package}
\alias{flowQB-package}
\alias{flowQB}
\docType{package}
\title{ Automatic Q and B Calculations}
\description{

We developed an automated approach to determine a cytometer's detection efficiency (Q) and background illumination (B) independent of a pre-defined bead set, based on Kmeans clustering and quadratic regression methods. We used a quadratic formulation for modelling Q and B that offers more physical insight about cytometry sensitivity than its truncated linear version, by considering the term CVintrinsic as part of the problem to solve. Using Kmeans in place of manual gating enabled an automated analysis to calculate  Q, B and CVintrinsic quantities objectively and in a time-efficient manner. This flowQB package is a R implementation of the theoretical description of Q, B and CVintrinsic calculations in our manuscript "Automated Quadratic Characterization of Flow Cytometer Instrument Sensitivity".

}

\details{
\tabular{ll}{
Package: \tab flowQB\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-11-15\cr
License: \tab  \cr
LazyLoad: \tab yes\cr
}

A number of generic functions are used to conduct an automatic Q, B and CVintrinsic calculations, the following functions are considered as introductory:


Function ReadDD reads the given FCS file and removes doublet events in the channel of interest(Chang). The ReadDD function returns a 2D array having the mean flourescent intensities(MFI) of the singlet events, first column is for the channel of interest and the second column  is for the  companion channel to be used to facilitate the identification of the bead sub-populations.

Function KmeansMeanSD takes the 2D MFI array and generates a number of clusters and return their MFI means and SDs. 

Function MFI2MESF converts the obtained MFI means to MESF means. For instance, the MFI output of KmeansMeanSD are converted to MESF values with an option to correct the SDs.

Functions lrMESF and qrMESF are the linear and quadratic regressions which use the obtained Means and Variances to calculate the Q and B values. 

Function DiscriminantExamination uses the values in the output of the function qrMESF to estimate the discriminant of the resulting quadratic equation and can be used as an additional interpretation tool to aid in understanding cytometer sensitivity.

An advanced processing is in the vignettes AdvancedflowQBNIH2 and AdvancedflowQBNIH3 and have a detailed illustration how to use the functions:  gPS, rPS, MultilevelgPS, MultilevelrPS and qrWEIGHTEDMESF. For this introductory processing, examples are provided here and in the vignette IntroductoryflowQBNIH.

}
\author{
Faysal El Khettabi
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

A. Gaigalas and L. Wang, Approaches to Quantitation in Flow Cytometry, in Standardization and Quality Assurance in Fluorescence Measurements II Springer Series on Fluorescence (2008), Volume 6, Part D, 371-398.

}

\examples{

rm(list=ls(all=TRUE))
library("flowQB")
File= system.file("extdata","NIH.fcs",package="flowQB")

# Doublet Discriminations 
  	DD=96
# Reading the FCS file, the Forward Scattering Area 1 and 2 
# are used to obtain singlet events with DD=96%.
# The processing returns a 2D singlet events for the channel 
# of interest 5 with the companion channel Side Scattering 3.
  	t2D=ReadDD(File,1,2,DD,5,3)
# Function KmeansMeanSD returns the MFI Means and SDs of the 8 clusters.
  	MeansSDs=KmeansMeanSD(t2D,8,500,200,1)
# CVs calculation
  	CVs=MeansSDs[,2]/MeansSDs[,1]
# MESF calculation, 
  	p=357217.00/7102
# Linear Q and B Calculation.
# MFIs are converted to MESFs and SDs are corrected using the beads in cluster 8.
  	l2D=MFI2MESF(MeansSDs,p,CVs[8])
# Linear regression 
# Peaks associated to cluster  3 to  cluster 5 are used 
# to compute the linear regression coefficients.
  	LQB=lrMESF(l2D,3,5)
  		print("Linear QB")
  			print(LQB)
# Quadratic Q and B Calculation.
# MFIs are converted to MESFs and SDs are NOT corrected.
  	Q2D=MFI2MESF(MeansSDs,p,0)
# Peaks associated to cluster  3 to  cluster 6 are used to compute 
# the quadratic regression coefficients.
  	QQB=qrMESF(Q2D,3,6)
  		print("Quadratic QB")
  			print(QQB)

}
