\name{fit_led}
\alias{fit_led}

\title{ 
    TODO
}

\description{
    TODO.
}

\usage{
    fit_led(fcs_file_path_list, ignore_channels, dyes, detectors, 
        bounds, signal_type, instrument_name, minimum_rows = 3, 
        max_iterations = 10, ...)
}

\arguments{
    \item{fcs_file_path_list}{
        TODO.
    }
    \item{ignore_channels}{
        TODO.
    }
    \item{dyes}{
        TODO.
    }
    \item{detectors}{
        TODO.
    }
    \item{bounds}{
        TODO.
    }
    \item{signal_type}{
        TODO.
    }
    \item{instrument_name}{
        TODO.
    }
    \item{minimum_rows}{
        TODO.
    }
    \item{max_iterations}{
        TODO.
    }
    \item{\ldots}{
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
    \code{\link{fit_multipeak}}
}

\examples{
    library(flowCore)
    library(xlsx)

    inst_xlsx_path <- system.file("extdata", "example1", 
        "140126_InstEval_Stanford_LSRIIA2.xlsx", package="flowQB")
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
    folder_column <- 18
    folder_row <- 14
    folder <- xlsx[[folder_column]][[folder_row]]
    fcs_directory <- system.file("extdata", "example1",
        instrument_folder, folder, package="flowQB")
    fcs_file_path_list <- list.files(fcs_directory, "*.fcs", full.names= TRUE)

    bounds_min_col <- 6
    bounds_min_row <- 7
    bounds_max_col <- 7
    bounds_max_row <- 7
    bounds <- list()
    if (is.na(xlsx[[bounds_min_col]][[bounds_min_row]])) {
        bounds$minimum <- -100
    } else {
        bounds$minimum <- as.numeric(xlsx[[bounds_min_col]][[bounds_min_row]])
    }
    if (is.na(xlsx[[bounds_max_col]][[bounds_max_row]])) {
        bounds$maximum <- 100000
    } else {
        bounds$maximum <- as.numeric(xlsx[[bounds_max_col]][[bounds_max_row]])
    }
    signal_type_col <- 3
    signal_type_row <- 19
    signal_type <- xlsx[[signal_type_col]][[signal_type_row]]
    
    instrument_name_col <- 2
    instrument_name_row <- 5
    instrument_name <- xlsx[[instrument_name_col]][[instrument_name_row]]

    channel_cols <- 3:12
    dye_row <- 11
    detector_row <- 13
    dyes <- as.character(xlsx[dye_row,channel_cols])
    detectors <- as.character(xlsx[detector_row,channel_cols])

    led_results <- fit_led(fcs_file_path_list, ignore_channels, dyes,
        detectors, bounds, signal_type, instrument_name, 
        minimum_rows = 3, max_iterations = 10)
}

\keyword{functions}
