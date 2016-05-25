\name{LEDflowQBCalculation}
\alias{LEDflowQBCalculation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QB Calculation for LEDs
}
\description{
This function is used for the LED FCS files to determine, for the channels of interest, the raw statitics for the regression 
and the generation of the regression's coefficients, Q and B values. 
}
\usage{
LEDflowQBCalculation(input, CHANNELTOBEPROCESSED, minF, maxF)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{input}{
Folder where the LED FCS files is stored.
}
  \item{CHANNELTOBEPROCESSED}{
Indices of channels of interest in the LED FCS files.
}
  \item{minF}{
Minimum MFI to consider for the peaks.
}
  \item{maxF}{
Maximum MFI to consider for the peaks.
}
}
\details{
The function generates the results as a list, the first element of the list is for raw statistics and the second element of the list is for the coefficients. 
}
\value{
Q and B values with their "Pvalue" & "Std-Error" are in the second list with the coefficients.
}
\references{ 
Faysal El Khettabi et al. 2014, Automated Quadratic Characterization of Flow Cytometer Instrument Sensitivity, to be submitted. }
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
