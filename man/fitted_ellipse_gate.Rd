\name{fitted_ellipse_gate}
\alias{fitted_ellipse_gate}

\title{ Singlets Detection}
\description{
Gate by an ellipse of radius R on the scatter data.
}
\usage{
fitted_ellipse_gate(myFlowFrame, channels, R=1)
}

\arguments{
  \item{myFlowFrame}{
FCS file.
}
  \item{channels}{
Indices of the Forward scatter and Side Scatter in the FCS file.
}
  \item{R}{
Gate by an ellipse of radius R on the scatter data.

}
}
\details{
Singlets detection using Forward scatter and Side Scatter in the FCS file.
}
\value{
A new FCS file with only singlets.
}

\author{
    Wayne Moore, Faysal El Khettabi, Josef Spidlen
}

\seealso{
    \code{\link{find_peak}}
}

\examples{
    library('flowCore')
    fcsFilePath <- system.file("extdata", "935289.fcs", package="flowQB")
    myFlowFrame <- read.FCS(fcsFilePath)
    gatedFlowFrame <- fitted_ellipse_gate(myFlowFrame, c('FSC-H', 'SSC-H'))
}

\keyword{functions}
\keyword{helper functions}

