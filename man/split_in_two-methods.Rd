\name{split_in_two}
\alias{split_in_two}
\alias{split_in_two,matrix,split_in_two-method}
\alias{split_in_two,matrix,ANY-method}
\alias{split_in_two,matrix-method}
\alias{split_in_two,flowFrame,split_in_two-method}
\alias{split_in_two,flowFrame,ANY-method}
\alias{split_in_two,flowFrame-method}

\title{Split an object in a low density region}

\description{
    This method finds a split in the low density region of the provided object
    (which shall be either a matrix or flowCore's flowFrame object) and 
    returns a vector of TRUE/FALSE depending on whether each of the events
    (rows) are left or righ (i.e, lower or higher) than 
    density peak. An FCS channel shall be specified if a flowFrame object
    with multiple channels in provided on the input. If matrix is the input
    then it shall contain one column only. This method is designed to work
    well for 2 level beads, such as stained and unstained, but will not
    provide meaningful results for multi-level beads or other data in
    general.
}

\usage{
    split_in_two(object, \dots)
}

\arguments{
    \item{object}{Object of class \code{\linkS4class{flowFrame}} or a matrix.}
    \item{\dots}{Additional options, see the details section.}
}

\details{
    Additional parameter of the method:
    \describe{
        \item{\code{channel}}{Which FCS channel shall be used in order to
        look for the split in the data? This is applicable if a flowFrame
        object with several channels is used on the input.}
    }
}

\value{
    A vector of TRUE/FALSE values depending on whether each of the events 
    (rows) are below or above the identified split value.
}

\author{
    Wayne Moore, Faysal El Khettabi, Josef Spidlen
}

\examples{
    library('flowCore')
    ## library(flowQBData)
    ## flowQBData is available since BioConductor 3.4, please install it
    ## manually in order to be able to these examples

#   fcsFilePath <- system.file("extdata", "SSFF_LSRII", "Other_Tests",
#       "933745.fcs", package="flowQBData")
#   myFlowFrame <- read.FCS(fcsFilePath)
#   ## Note that this is just to demonstrate the syntax, but doing this on the
#   ## FSC and SSC channels of this particular FCS file is not very meaningful
#   r1 <- split_in_two(myFlowFrame, 'FSC-H')
#   r2 <- split_in_two(exprs(myFlowFrame[,'SSC-H']))
}

\keyword{methods}
