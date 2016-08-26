\name{calc_mean_sd_background}
\alias{calc_mean_sd_background}

\title{ 
    Calculate the mean and the standard deviation.
}

\description{
    This method looks at all channels except those specified in the
    ignore_channel list and calculated the mean and the standard deviation for
    those channels.
}

\usage{
    calc_mean_sd_background(fcs_file_path, ignore_channels)
}

\arguments{
    \item{fcs_file_path}{
A character string specifying the file path
to the FCS file with the acquired bead data.    
    }
    \item{ignore_channels}{
A vector of short channel names (values of
the $PnN keywords) specifying channels that should not be considered for
the fitting procedure. Normally, those should be all non-fluorescence channels,
such as the time and the (forward and side) scatter channels.
    }
}

\details{
    The getOutliers method from the extremevalues package is used to
    calculate the mean and the standard deviation values for
    all the FCS channels in the file except those specified in the
    ignore channels list.
}

\value{
    The result is a data frame with columns corresponding to short channel
    names of channels from the input FCS file except those specified by the
    ignore_channels parameter. The rows include the total number of events,
    the mean and the standard deviation.
}

\author{
    Josef Spidlen, Wayne Moore, Faysal El Khettabi
}

\seealso{
    \code{\link{calc_mean_sd_197}},
    \code{\link{calc_mean_sd_duke}},
    \code{\link{calc_mean_sd_capture}},
    \code{\link{calc_mean_sd_capture_all}}
}

\examples{
    library(flowCore)
    library(xlsx)
    library(flowQBData)

    inst_xlsx_path <- system.file("extdata", 
        "140126_InstEval_Stanford_LSRIIA2.xlsx", package="flowQBData")
    xlsx <- read.xlsx(inst_xlsx_path, 1, headers=FALSE, stringsAsFactors=FALSE)
    
    ignore_channels_row <- 9
    ignore_channels <- vector()
    i <- 1
    while(!is.na(xlsx[[i+4]][[ignore_channels_row]])) {
        ignore_channels[[i]] <- xlsx[[i+4]][[ignore_channels_row]]
        i <- i + 1
    }
    
    instrument_folder_row <- 9
    instrument_folder_col <- 2
    instrument_folder <- xlsx[[instrument_folder_col]][[instrument_folder_row]]

    test_column <- 15
    test_row <- 14
    folder <- xlsx[[test_column]][[test_row]]
    file_name <- xlsx[[test_column]][[test_row+1]]

    fcs_path <- system.file("extdata",
        instrument_folder, folder, file_name, package="flowQBData")

    results <- calc_mean_sd_background(fcs_path, ignore_channels)

    ## Same thing as above with providing the arguments directly without
    ## parsing it from the spreadsheet.
    fcs_path <- system.file("extdata", "SSFF_LSRII", "Other_Tests", 
        "935319.fcs", package="flowQBData")
    ignore_channels <- c(
        "Time", "FSC-A", "FSC-W", "FSC-H", "SSC-A", "SSC-W", "SSC-H")
    results <- calc_mean_sd_background(fcs_path, ignore_channels)
}

\keyword{functions}
