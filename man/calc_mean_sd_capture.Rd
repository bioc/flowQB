\name{calc_mean_sd_capture}
\alias{calc_mean_sd_capture}

\title{ 
    TODO
}

\description{
    TODO.
}

\usage{
    calc_mean_sd_capture(fcs_file_path, scatter_channels, detector, dye)
}

\arguments{
    \item{fcs_file_path}{
        TODO.
    }
    \item{scatter_channels}{
        TODO.
    }
    \item{detector}{
        TODO.
    }
    \item{dye}{
        TODO.
    }
}

\details{
    TODO.
}

\value{
    TODO.
}

\author{
    Wayne Moore, Faysal El Khettabi, Josef Spidlen
}

\seealso{
    \code{\link{calc_mean_sd_capture_all}}
}

\examples{
    library(flowCore)
    library(flowQBData)

    fcs_file_path <- system.file("extdata", "SSFF_LSRII", "SU_2B",
        "933723.fcs", package="flowQBData")

    scatter_channels <- c("FSC-A", "SSC-A")
    detector <- "APC-A"
    dye <- "APC"
    
    results <- calc_mean_sd_capture(
        fcs_file_path, scatter_channels, detector, dye)
}

\keyword{functions}

