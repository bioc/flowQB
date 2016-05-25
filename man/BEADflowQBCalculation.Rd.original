\name{BEADflowQBCalculation}
\alias{BEADflowQBCalculation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QB Calculation for Beads
}
\description{
This function is used for the bead FCS file to determine the singlet events. These singlet events are clustered for the channels of interest to determine the raw statitics for the regression and the generation of the regression's coefficients, Q and B values.  
}
\usage{
BEADflowQBCalculation(file, R, width, fraction, numberoftotalepeaks, scatters, CHANNELTOBEPROCESSED, pi, pf, Viz,output)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{file}{
Name of the FCS file.
}
  \item{R}{
Gate by an ellipse of radius R on the scatter data.
}
  \item{width}{
Parameter used to find the lower bound (default width at half maximum) and the upper bound (default .5).
}
  \item{fraction}{
Parameter used to find the densest fraction of cells (default .1).
}
  \item{numberoftotalepeaks}{
Global Number of the peaks in the data.
}
  \item{scatters}{
Indices of Forward Scatter (FS) and Side scatter (SS) in the bead FCS file.
}
  \item{CHANNELTOBEPROCESSED}{
Indices of channels of interest in the bead FCS file.
}
  \item{pi}{
Dimmest peak to be used in the regression.
}
  \item{pf}{
Last peak to be used in the regression.
}
  \item{Viz}{
 Parameter to control the visualization of Kmeans results only if Viz=1, a subfolder called flowQB in the output folder has subfolder "PLOTSKMEANS_BEADS_" with the visualizations.
}
  \item{output}{ 
Folder where the Kmeans image results should be stored.
}
}
\details{
The function generates the results as a list, the first is for RAW STATISTICS and the second is for THE COEFFICIENTS. 
}
\value{
Q and B values with their "Pvalue" & "Std-Error" are in the second list with the coefficients.
}
\references{ 
Faysal El Khettabi et al. 2014, Automated Quadratic Characterization of Flow Cytometer Instrument Sensitivity, to be submitted.
}
\author{
Faysal El Khettabi et al.
}

\seealso{ 

J. Wood, Fundamental Flow Cytometer Properties Governing Sensitivity and Resolution, 
Cytometry 33, (1998), p.~ 260 - 6.

E. Chase and R. Hoffman, Resolution of Dimly Fluorescent Particles: a Practical Measure of Fluorescence Sensitivity,
Cytometry 33 (1998), p.~ 267-279.

R. Hoffman and J. Wood, Characterization of Flow Cytometer Instrument Sensitivity, Current Protocols in Cytometry, Chapter 1: Unit 1.20 (2007).

A. Gaigalas and L. Wang, Approaches to Quantitation in Flow Cytometry, in Standardization and Quality Assurance in Fluorescence Measurements II Springer Series on Fluorescence (2008), Volume 6, Part D, 371-398.

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
