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
    
    ## Now the same thing again, but we will show how to extract information
    ## from the spreadsheet and run the appropriate calculations
    library(xlsx)
    xls_path <- system.file("extdata", "example1", 
        "140126_InstEval_Stanford_LSRIIA2.xlsx", package="flowQB")
    xls <- read.xlsx(xls_path, 1, headers=FALSE, stringsAsFactors=FALSE)
    insfolder <- instrument.folder <- xls[[2]][[9]]

    dyes <- list()
    detectors <- list()
    filepaths <- list()
    scatters <- list()

    for (i in 1:10)
    {
        folder <- xls[[i+2]][[14]]
        filename <- xls[[i+2]][[15]]

        if (is.na(filename)) next
        filepath <- system.file("extdata", "example1", insfolder, 
            folder, filename, package="flowQB")
        ## Spreadsheet may describe additional FCS files not included
        ## with the library, so skip if file doesn't exist
        if (nchar(filepath) == 0) next

        filepaths <- c(filepaths, filepath)
        dyes <- c(dyes, xls[[i+2]][[11]])
        detectors <- c(detectors, xls[[i+2]][[13]])
        scatters[[length(scatters)+1]] <- c(xls[[i+2]][[16]], xls[[i+2]][[17]])
    }

    results2 <- calc_mean_sd_capture_all(filepaths, scatters, detectors, dyes)
}

\keyword{functions}

