\name{calc_mean_sd_capture_all}
\alias{calc_mean_sd_capture_all}

\title{ 
    TODO
}

\description{
    TODO.
}

\usage{
    calc_mean_sd_capture_all(fcs_file_path_list, scatter_channels_list,
    detector_list, dye_list)
}

\arguments{
    \item{fcs_file_path_list}{
        TODO.
    }
    \item{scatter_channels_list}{
        TODO.
    }
    \item{detector_list}{
        TODO.
    }
    \item{dye_list}{
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
    Josef Spidlen, Wayne Moore, Faysal El Khettabi
}

\seealso{
    \code{\link{calc_mean_sd_capture}}
}

\examples{
    library(flowCore)

    file_directory <- system.file("extdata", "example1", 
        "SSFF_LSRII", "SU_2B", package="flowQB")
    fcs_file_path_list <- as.list(file.path(
        file_directory, c("933723.fcs","933725.fcs")))
    scatter_channels_list <- list(c("FSC-A", "SSC-A"), c("FSC-A", "SSC-A"))
    detector_list <- list("APC-A", "APC-Cy7-A")
    dye_list <- list("APC", "APC-Cy7")
    
    results <- calc_mean_sd_capture_all(
        fcs_file_path_list, 
        scatter_channels_list, 
        detector_list, 
        dye_list
    )
}

\keyword{functions}

