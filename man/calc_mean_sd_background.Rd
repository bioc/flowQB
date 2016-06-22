\name{calc_mean_sd_background}
\alias{calc_mean_sd_background}

\title{ 
    TODO
}

\description{
    TODO.
}

\usage{
    calc_mean_sd_background(fcs_file_path, ignore_channels)
}

\arguments{
    \item{fcs_file_path}{
        TODO.
    }
    \item{ignore_channels}{
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
