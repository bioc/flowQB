\name{pick_parameters}
\alias{pick_parameters}
\alias{pick_parameters,flowFrame,pick_parameters-method}
\alias{pick_parameters,flowFrame,ANY-method}
\alias{pick_parameters,flowFrame-method}

\title{Pick channel names}

\description{
Extract all channel names from a flowFrame object or column names from 
a matrix except those specified in a provided ignore list.
}

\usage{
    pick_parameters(object, ignore)
}

\arguments{
    \item{object}{Object of class \code{\linkS4class{flowFrame}} or a matrix.}
    \item{ignore}{A vector of channel names that we want to ignore}
}

\details{
    This method simply looks at all channel names in the provided
    flowFrame object or all column names of a matrix, then subtracts those
    specified in the ignore list and returns the resulting vector of
    channel/column names.
}

\value{
    A vector of character strings containing channel names of channels
    that were in the input flowFrame object (or columns of the input matrix)
    but were not included in the provided ignore list.
}

\author{
    Josef Spidlen, Wayne Moore, Faysal El Khettabi
}

\examples{
    library('flowCore')
    ## library(flowQBData)
    ## flowQBData is available since BioConductor 3.4, please install it
    ## manually in order to be able to these examples

#   fcsFilePath <- system.file("extdata", "SSFF_LSRII", "Other_Tests",
#       "933745.fcs", package="flowQBData")
#   myFlowFrame <- read.FCS(fcsFilePath)
#   ignore <- c("Time", "FSC-H", "FSC-A", "FSC-W", "SSC-H", "SSC-A", "SSC-W")
#   fluorescences <- pick_parameters(myFlowFrame, ignore)
}

\keyword{methods}
