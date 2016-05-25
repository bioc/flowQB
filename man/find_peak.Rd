\name{find_peak}
\alias{find_peak}

\title{ 
    Peak Detection
}

\description{
    The determination of the peak uses the densest fraction of events,
    the lower bound and the upper bound.
}
\usage{
    find_peak(data, width = 0.5, fraction = 0.1)
}

\arguments{
    \item{data}{
        FCS file.
    }
    \item{width}{
        Parameter used to find the lower bound (default width at half maximum)
        and the upper bound (default .5).
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
    Wayne Moore et al. 2014, Automated Quadratic Characterization of Flow
    Cytometer Instrument Sensitivity, to be submitted.
}

\author{
    Wayne Moore
}

\examples{
    my_peak_info <- find_peak(rnorm(1000, mean=5))
}

\keyword{functions}
\keyword{helper functions}
