\name{calc_mean_sd_capture_all}
\alias{calc_mean_sd_capture_all}

\title{ 
    Calculate the mean and the standard deviation 
    for the stained and unstained population
    of specified channels in specified FCS files.
}

\description{
    This methods performs the \code{calc_mean_sd_capture} function
    on a list of FCS files, list of scatter channel pairs, list of detectors
    and a list of dyes, and collates the results. The order of the arguments
    in the input lists matters, i.e., the first FCS file will be matched
    with the first pair of FSC/SSC channel names, the first detector
    name and the first dye name.
}

\usage{
    calc_mean_sd_capture_all(fcs_file_path_list, scatter_channels_list,
    detector_list, dye_list)
}

\arguments{
    \item{fcs_file_path_list}{
        A list of n FCS files, one for each detector.
    }
    \item{scatter_channels_list}{
        A list of n pairs of forward and side scatter channel names.
    }
    \item{detector_list}{
        A list of n detector names; those shall correspond to specific
        detector in the n specified FCS files.
    }
    \item{dye_list}{
        A list of n dye names; those will be used to name the columns of the
        resulting data frame.
    }
}

\details{
    This method assumes that each of the FCS files have useful data 
    only in the specified channel. Therefore, we perform the 
    \code{calc_mean_sd_capture} on all these FCS files separatelly and
    then put the results together into a single data frame.
}

\value{
    The result is a data frame with n columns, the headings of the
    columns correspond to the values in the list provided by the dye_list 
    argument. The rows include the
    total number of events,  the number of events in the FSC/SSC ellipse gate,
    the number of events in the high peak gate and low peak gate, the stained
    mean and stained standard deviation (based on the high peak gate), and
    finally the unstained mean and unstained standard deviation (based on
    the low peak gate).
}

\author{
    Josef Spidlen, Wayne Moore, Faysal El Khettabi
}

\seealso{
    \code{\link{calc_mean_sd_capture}}
}

\examples{
    library(flowCore)
    library(flowQBData)

    file_directory <- system.file("extdata", "SSFF_LSRII", "SU_2B", 
        package="flowQBData")
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
    xls_path <- system.file("extdata", "140126_InstEval_Stanford_LSRIIA2.xlsx",
        package="flowQBData")
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
        filepath <- system.file("extdata", insfolder, folder, filename, 
            package="flowQBData")
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

