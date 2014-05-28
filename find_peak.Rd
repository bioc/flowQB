\name{find_peak}
\alias{find_peak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Peak Detection
}
\description{
The determination of the peak uses the densest fraction (default fraction = .1 ) of events, the lower bound and the upper bound (default width = .5 )
}
\usage{
find_peak(data, width = 0.5, fraction = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
FCS file.
}
  \item{width}{
Parameter used to find the lower bound (default width at half maximum) and the upper bound (default .5).
}
  \item{fraction}{
Parameter used to find the densest fraction of cells (default .1).
}
}
\details{
This function will be used with the function scatter_gate.
}
\value{
Peaks for the scatters.
}
\references{
Wayne Moore et al. 2014, Automated Quadratic Characterization of Flow Cytometer Instrument Sensitivity, to be submitted. 
}
\author{
Wayne Moore
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
