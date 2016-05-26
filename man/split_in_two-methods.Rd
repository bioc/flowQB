\name{split_in_two}
\alias{split_in_two}
\alias{split_in_two,matrix,split_in_two-method}
\alias{split_in_two,matrix,ANY-method}
\alias{split_in_two,matrix-method}
\alias{split_in_two,flowFrame,split_in_two-method}
\alias{split_in_two,flowFrame,ANY-method}
\alias{split_in_two,flowFrame-method}

\title{TODO}

\description{
    TODO
}

\usage{
    split_in_two(object, \dots)
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
    fcsFilePath <- system.file("extdata", "935289.fcs", package="flowQB")
    myFlowFrame <- read.FCS(fcsFilePath)
    r1 <- split_in_two(myFlowFrame, 'FSC-H')
    r2 <- split_in_two(exprs(myFlowFrame[,'SSC-H']))
}

\keyword{methods}
