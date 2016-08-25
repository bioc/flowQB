\name{peak_gate}
\alias{peak_gate}
\alias{peak_gate,matrix,peak_gate-method}
\alias{peak_gate,matrix,ANY-method}
\alias{peak_gate,matrix-method}
\alias{peak_gate,flowFrame,peak_gate-method}
\alias{peak_gate,flowFrame,ANY-method}
\alias{peak_gate,flowFrame-method}

\title{Gate a 1D density peak in the provided object}

\description{
    This method finds a density peak in the provided object (which shall be
    either a matrix or flowCore's flowFrame object) and returns a vector of
    TRUE/FALSE depending on whether each of the events (rows) are in the
    density peak. An FCS channel shall be specified if a flowFrame object
    with multiple channels in provided on the input.
}

\usage{
    peak_gate(object, \dots)
}

\arguments{
    \item{object}{Object of class \code{\linkS4class{flowFrame}} or a matrix.}
    \item{\dots}{Additional options, see the details section.}
}

\details{
    TODO

    \describe{
        \item{\code{channel}}{Which FCS channel shall be used in order to
        look for the density peak? This is applicable if a flowFrame
        object with several channels is used on the input.}
        \item{\code{R}}{The radius to be used when finding the peak; R=1 by
        default.}
    }
}

\value{
    A vector of TRUE/FALSE values depending on whether each of the events 
    (rows) are located in the identified density peak.
}

\author{
    Wayne Moore, Faysal El Khettabi, Josef Spidlen
}

\examples{
    library('flowCore')
    library('flowQBData')
    
    fcsFilePath <- system.file("extdata", "SSFF_LSRII", "Other_Tests",
        "933745.fcs", package="flowQBData")

    myFlowFrame <- read.FCS(fcsFilePath)
    r1 <- peak_gate(myFlowFrame, 'FSC-H')
    r2 <- peak_gate(exprs(myFlowFrame[,'SSC-H']))
    ## r3 will have more events than r2
    r3 <- peak_gate(exprs(myFlowFrame[,'SSC-H']), R=1.5)
}

\keyword{methods}
