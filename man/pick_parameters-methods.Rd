\name{pick_parameters}
\alias{pick_parameters}
\alias{pick_parameters,flowFrame,pick_parameters-method}
\alias{pick_parameters,flowFrame,ANY-method}
\alias{pick_parameters,flowFrame-method}

\title{TODO}

\description{
    TODO
}

\usage{
    pick_parameters(object, ignore)
}

\arguments{
    \item{object}{Object of class \code{\linkS4class{flowFrame}} or a matrix.}
    \item{ignore}{TODO.}
    %\item{\dots}{Additional options, see the details section.}
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
    Josef Spidlen, Wayne Moore, Faysal El Khettabi
}

\examples{
    library('flowCore')
    fcsFilePath <- system.file("extdata", "935289.fcs", package="flowQB")
    myFlowFrame <- read.FCS(fcsFilePath)
    ignore <- c("Time", "FSC-H", "FSC-A", "FSC-W", "SSC-H", "SSC-A", "SSC-W")
    fluorescences <- pick_parameters(myFlowFrame, ignore)
}

\keyword{methods}
