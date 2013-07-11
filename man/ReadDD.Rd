\name{ReadDD}
\alias{ReadDD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ReadDD
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function ReadDD reads a given FCS file and remove doublet events in the channel of interest(ChanGiven).
}
\usage{
ReadDD(File, Chan1DD, Chan2DD, P, ChanGiven, ChanCompanion)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{File}{
%%     ~~Describe \code{File} here~~
Path of the file to be read.
}
  \item{Chan1DD}{
%%     ~~Describe \code{Chan1DD} here~~
Number(index) of the first channel to be used in detecting singlet events.
}
  \item{Chan2DD}{
%%     ~~Describe \code{Chan2DD} here~~
Number(index) of the second channel to be used in detecting singlet events.
}
  \item{P}{
%%     ~~Describe \code{P} here~~
Doublet discrimination level, only events satistfying P/100 < Chan1DD/Chan2DD <= 2-P/100 are extracted by the function ReadDD and they are singlet events. P can be a number between 50 to 99.
}
 \item{ChanGiven}{
%%     ~~Describe \code{ChanGiven} here~~
Number(index) of the channel of interest.
}
  \item{ChanCompanion}{
%%     ~~Describe \code{ChanCompanion} here~~
Number(index) of the companion channel to be used to facilitate the identification of the bead sub-populations.
}

}

\value{
The ReadDD function returns a 2D array having the mean flourescent intensities(MFI) of the singlet events, first column is for the channel of interest(ChanGiven) and the second column  is for the  companion channel(ChanCompanion) to be used to facilitate the identification of the bead sub-populations.

}

\details{
The returned 2D array will be clustered using the function KmeansMeanSD.
}

\references{
%% ~put references to the literature/web site here ~
Faysal El Khettabi et al. 2013, Automated Quadratic Characterization of Flow Cytometer Instrument Sensitivity, to be submitted.
}
\author{Faysal El Khettabi \\
Terry Fox Laboratory \\
British Columbia Cancer Agency \\
Vancouver, BC, Canada \\
E-mail: fkhettabi@bccrc.ca 
}
\note{
It reads only a FCS file. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
flowCore package.
}
\examples{

rm(list=ls(all=TRUE))
rm(list=ls(all=TRUE))
library("flowQB")
File= system.file("extdata","NIH.fcs",package="flowQB")

#  Now FCS file is in data folder as a temp file.
#  This file will be removed at the end of the illustration.

# Doublet Discriminations 
P=96

# Reading the FCS file, the Forward Scattering Area 1 and 2 
# are used to obtain singlet events with P=96%.
# The processing returns a 2D singlet events for the channel 
# of interest, index 5, with the companion channel Side Scattering, index 3.
  
# Usage of the function ReadDD	
t2D=ReadDD(File,1,2,P,5,3)
summary(t2D)


}






