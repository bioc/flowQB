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
Wayne Moore
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
