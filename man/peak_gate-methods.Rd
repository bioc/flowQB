\name{peak_gate}
\alias{peak_gate}
\alias{peak_gate,matrix,peak_gate-method}
\alias{peak_gate,matrix,ANY-method}
\alias{peak_gate,matrix-method}
\alias{peak_gate,flowFrame,peak_gate-method}
\alias{peak_gate,flowFrame,ANY-method}
\alias{peak_gate,flowFrame-method}

\title{TODO}

\description{
    TODO
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
        \item{\code{channel}}{The path to the directory where to 
            download the files. By default, the files will be downloaded to a 
            folder named based on the dataset identifier (FR-FCM-xxxx) 
            that will be created in your working directory.}
    }
}

\value{
    TODO
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
}

\keyword{methods}
