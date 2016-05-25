\name{flowQB-package}
\alias{flowQB-package}
\alias{flowQB}
\docType{package}
\title{Automated Quadratic Characterization of Flow Cytometer Instrument Sensitivity: Q, B and CV-instrinsic calculations.
}
\description{
flowQB is a fully automated R Bioconductor package to calculate automatically the detector efficiency (Q), optical background (B) and intrinsic CV of the beads.
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

BEADflowQBCalculation: This function is used for the bead FCS file to determine the singlet events. These singlet events are clustered for the channels of interest to determine the raw statitics for the regression and the generation of the regression's coefficients, Q and B values.  

LEDflowQBCalculation: This function is used for the LED FCS files to determine, for the channels of interest, the raw statitics for the regression and the generation of the regression's coefficients, Q and B values. 

These functions generate the results as a list, the first element of the list is for Raw Statistics and the second element of the list is for the coefficients, Q and B values. 
 
Raw Statistics uses 3 approaches, Robust Statistics, Density estimation assuming a Gaussian distribution ( MASS package) and the second 'extremevalues' package used to determine the raw statistics without outliers: 
1) Number of events in each peak: NE  
2) MFI associated to Robust Statistics: mfiRS
3) MFI associated to Gauss estimation using MASS: mfiGS        
4) MFI associated to Gauss estimation using 'extremevalues': mfino
5) Standard deviation associated to Robust Statistics: mfirSD
6) Standard deviation associated to Gauss estimation using MASS: mfiGS 
6) Standard deviation associated to Gauss estimation using 'extremevalues': mfiSDno
7) Number of events in each peak without outliers: Nesno.ORD

For each approach (StatsProcedure) and channel (MARKER), the coefficients (c0,c1,c2) and (Q,B) are listed with their associated (Pvalue, Std-Error). 

}
\author{
Faysal EL Khettabi and Wayne Moore
}
\references{ 
Faysal El Khettabi et al. 2014, Automated Quadratic Characterization of Flow Cytometer Instrument Sensitivity, to be submitted. }
\keyword{ package }
\seealso{
J. Wood, Fundamental Flow Cytometer Properties Governing Sensitivity and Resolution, 
Cytometry 33, (1998), p.~ 260 - 6.
E. Chase and R. Hoffman, Resolution of Dimly Fluorescent Particles: a Practical Measure of Fluorescence Sensitivity,
Cytometry 33 (1998), p.~ 267-279.
R. Hoffman and J. Wood, Characterization of Flow Cytometer Instrument Sensitivity, Current Protocols in Cytometry, Chapter 1: Unit 1.20 (2007).
A. Gaigalas and L. Wang, Approaches to Quantitation in Flow Cytometry, in Standardization and Quality Assurance in Fluorescence Measurements II Springer Series on Fluorescence (2008), Volume 6, Part D, 371-398.
}
\examples{
if(1==0)
{
rm(list=ls(all=TRUE))
library(flowCore)
# TO READ FCS FILE, the flowcore library is needed
# Where is the file and its name? 
fpath = system.file("extdata","FCSSMALLEVENTS.fcs",package="flowQB")
# fpath="/home/fkhettabi/Desktop/newflowQB2014/flowQB/inst/extdata/FCSSMALLEVENTS.fcs"
# Name of the FCS file
file=fpath
# Singlet Parameter:
# R: gate by an ellipse of radius R on the scatter data
R=2
# width: a parameter used to find the lower bound (default width at half maximum) and the upper bound (default .5).
width=.5
# fraction: a parameter used to find the densest fraction of cells (default .1).
fraction=.1
# Number of the peaks in the data
numberoftotalepeaks=8
# To define the indices of the FS and SS.
# NOW GO THE THE INPUT 
fcs <- read.FCS(file)
names(fcs)
# Scatters FS & SS are:
scatters=c(2,5)
# Channels of interest, you need to select at least two channels, for instance:
CHANNELTOBEPROCESSED=c(8,9,36,37,38,39)
names(fcs)[CHANNELTOBEPROCESSED]
nClusters=numberoftotalepeaks
###################################################
### PEAKS TO BE USED IN THE REGRESSION
###################################################
pi=1
pf=7
output=""
RESULTS=BEADflowQBCalculation(file,R,width,fraction,numberoftotalepeaks,scatters,CHANNELTOBEPROCESSED,pi,pf,0,output)
}

if( 1==0)
{
rm(list=ls(all=TRUE))
NAMEOFTHEFILES=c("935295.fcs","935297.fcs","935299.fcs","935301.fcs","935303.fcs","935305.fcs")
# GET THE FILES WITH THEIR PATHS
# bulid "input"
input=rep("",length(NAMEOFTHEFILES))
for( i in 1:length(NAMEOFTHEFILES))
{
input[i] = system.file("extdata",NAMEOFTHEFILES[i],package="flowQB")
}

# TO READ FCS FILE, the flowcore library is needed
library(flowCore)
# read one file to have the name of the channels of interest!
fcs <- read.FCS(input[1])
orgdata=data.frame(exprs(fcs))
# Channels of interest, for instance:
CHANNELTOBEPROCESSEDx=c(8,9,36,37,38,39)
CHANNELTOBEPROCESSED = names(orgdata)[CHANNELTOBEPROCESSEDx]
# Minimum MFI to consider for the events.
minF=0
# Maximum MFI to consider for the events.
maxF=100000
RESULTS=LEDflowQBCalculation(input,CHANNELTOBEPROCESSED,minF,maxF)
}






}
