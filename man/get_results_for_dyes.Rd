\name{get_results_for_dyes}
\alias{get_results_for_dyes}

\title{ 
    Extract dye results from a data frame with detector results.
}

\description{
    This function takes a data frame where columns are named based on
    detectors and extracts a subset of the data frame it by selecting only
    specified detectors. In addition, the columns will be renamed based on the
    specified dyes argument.
}

\usage{
    get_results_for_dyes(dyes, detectors, results)
}

\arguments{
    \item{dyes}{
        A vector of n dye names which shall correspond to the dyes
        specified in the dyes argument. These will be the column names
        of the resulting data frame. The detector-dye mapping is done
        based on the order of values in the two vectors, i.e., the
        first dye shall correspond to the first detector, etc.
    }
    \item{detectors}{
        A vector of n detector names which shall correspond to the dyes
        specified in the dyes argument. These shall correspond to the
        column names in the input data frame. The detector-dye mapping is done
        based on the order of values in the two vectors, i.e., the
        first dye shall correspond to the first detector, etc.
    }
    \item{results}{
        An input data frame that shall contain columns corresponding to all
        the different values specified by the detectors vector.
    }
}

\details{
    This function is used to select a subset of columns from a data frame
    by specifying the columns of interest (detectors). In addition, the
    columns will be renamed to dyes corresponding to those detectors.
}

\value{
    A data frame with n columns, column names corresponding to the specified
    dyes and rows/values extracted from the input data frame.
}

\author{
    Wayne Moore, Faysal El Khettabi, Josef Spidlen
}

\seealso{
    \code{\link{calc_mean_sd_duke}}
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

    test_column <- 13
    test_row <- 14
    folder <- xlsx[[test_column]][[test_row]]
    beads_file_name <- xlsx[[test_column]][[test_row+1]]
    scatter_channels <- c(
        xlsx[[test_column]][[test_row+2]], 
        xlsx[[test_column]][[test_row+3]])

    fcs_path <- system.file("extdata",
        instrument_folder, folder, beads_file_name, package="flowQBData")

    results <- calc_mean_sd_duke(fcs_path, scatter_channels, ignore_channels)
    
    channel_cols <- 3:12
    dye_row <- 11
    detector_row <- 13
    dyes <- as.character(xlsx[dye_row,channel_cols])
    detectors <- as.character(xlsx[detector_row,channel_cols])
    dye_results <- get_results_for_dyes(dyes, detectors, results)
}

\keyword{functions}
\keyword{helper functions}
