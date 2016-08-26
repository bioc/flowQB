\name{fitted_ellipse_gate}
\alias{fitted_ellipse_gate}
\alias{fitted_ellipse_gate,flowFrame,fitted_ellipse_gate-method}
\alias{fitted_ellipse_gate,flowFrame,ANY-method}
\alias{fitted_ellipse_gate,flowFrame-method}

\title{Fit and apply an ellipse gate}
\description{
    Fit an ellipse (or ellipsoid) gate on the most dense region of the selected
    channels of a flowFrame object.
}

\usage{
    fitted_ellipse_gate(object, channels, R=1)
}

\arguments{
    \item{object}{
        A flowFrame object that will be gated by an automatically fitted
        ellipse (or ellipsoid) gate.
    }
    \item{channels}{
        A vector if indices or short channel names (valued of the 
        \$PnN keywords) of channels that shall be used to fit an ellipse
        (or ellipsoid) gate. Typically, these would be the indeces or names
        of the forward and scatter channels in the input flowFrame, which will
        then be gated by a 2D ellipse gate (since 2 channels were provided).
        A multidimensional ellipsoid gate is fitted analogically if more than
        2 channels are provided. A 1D "ellipsoid" gate (i.e., a range gate)
        is applied if only a single channel is specified.
    }
    \item{R}{
        An additional scaling factor for the radii of the ellipse/ellipsoid 
        gate.
    }
}

\details{
    First, events with negative values in the specified channels are removed.
    Next, the \code{\link{find_peak}} function is called on log transformed
    values of each of the specified channels. The returned lower and upper
    boundaries of each of the channels are used to fit an ellipse/ellipsoid
    gate. The radii of this gate are scalled by the value of the provided
    R argument. A new flowFrame object with only those events in this gate
    is returned.
}

\value{
    A flowFrame object gated by an automatically fitted ellipse or ellipsoid 
    gate created in the specified channels.
}

\author{
    Wayne Moore, Faysal El Khettabi, Josef Spidlen
}

\seealso{
    \code{\link{find_peak}}
}

\examples{
    library('flowCore')
    library('flowQBData')
    
    fcsFilePath <- system.file("extdata", "SSFF_LSRII", "Other_Tests",
        "933745.fcs", package="flowQBData")
    myFlowFrame <- read.FCS(fcsFilePath)
    gatedFlowFrame <- fitted_ellipse_gate(myFlowFrame, c('FSC-H', 'SSC-H'))
}

\keyword{functions}
\keyword{helper functions}
