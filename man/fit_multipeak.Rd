\name{fit_multipeak}
\alias{fit_multipeak}
\alias{fit_spherotech}
\alias{fit_thermo_fischer}

\title{ 
    TODO
}

\description{
    TODO.
}

\usage{
    fit_multipeak(fcs_file_path, scatter_channels, ignore_channels,
        N.peaks, dyes, detectors, bounds,
        signal_type, instrument_name, minimum_useful_peaks = 3, 
        max_iterations = 10, logicle_width = 0.5, ...)
    
    fit_spherotech(fcs_file_path, scatter_channels, ignore_channels,
        dyes, detectors, bounds,
        signal_type, instrument_name, minimum_useful_peaks = 3, 
        max_iterations = 10, logicle_width = 0.5, ...)

    fit_thermo_fischer(fcs_file_path, scatter_channels, ignore_channels,
        dyes, detectors, bounds,
        signal_type, instrument_name, minimum_useful_peaks = 3, 
        max_iterations = 10, logicle_width = 0.5, ...)
}

\arguments{
    \item{fcs_file_path}{
        TODO.
    }
    \item{scatter_channels}{
        TODO.
    }
    \item{ignore_channels}{
        TODO.
    }
    \item{N.peaks}{
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
    \item{minimum_useful_peaks}{
        TODO.
    }
    \item{max_iterations}{
        TODO.
    }
    \item{logicle_width}{
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
    \code{\link{fit_led}}
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
    
    folder_column <- 16
    folder_row <- 14
    folder <- xlsx[[folder_column]][[folder_row]]
    filename <- xlsx[[folder_column]][[folder_row+1]]
    scatter_channels <- c(
        xlsx[[folder_column]][[folder_row+2]], 
        xlsx[[folder_column]][[folder_row+3]])
        
    fcs_file_path <- system.file("extdata", instrument_folder, folder, 
        filename, package="flowQBData")

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

    multipeak_results <- fit_spherotech(fcs_file_path, scatter_channels, 
        ignore_channels, dyes, detectors, bounds, 
        signal_type, instrument_name, minimum_useful_peaks = 3, 
        max_iterations = 10, logicle_width = 0.5)

    ## The above is the same as this:
    ## N.peaks <- 8
    ## multipeak_results <- fit_multipeak(fcs_file_path, scatter_channels, 
    ##     ignore_channels, N.peaks, dyes, detectors, bounds, 
    ##     signal_type, instrument_name, minimum_useful_peaks = 3, 
    ##     max_iterations = 10, logicle_width = 0.5)

    plot(
        exprs(multipeak_results$transformed_data[,"FITC-A"]),
        exprs(multipeak_results$transformed_data[,"Pacific Blue-A"]),
        col=multipeak_results$peak_clusters$cluster, pch='.')
        
    ## Thermo-Fisher Example:
    folder_column <- 17
    folder <- xlsx[[folder_column]][[folder_row]]
    filename <- xlsx[[folder_column]][[folder_row+1]]

    fcs_file_path <- system.file("extdata", instrument_folder, folder, 
        filename, package="flowQBData")

    multipeak_results_tf <- fit_thermo_fischer(fcs_file_path, scatter_channels,
        ignore_channels, dyes, detectors, bounds, 
        signal_type, instrument_name, minimum_useful_peaks = 3, 
        max_iterations = 10, logicle_width = 0.5)
        
    ## The above is the same as this:
    ## N.peaks <- 6
    ## multipeak_results_tf <- fit_multipeak(fcs_file_path, scatter_channels,
    ##     ignore_channels, N.peaks, dyes, detectors, bounds, 
    ##     signal_type, instrument_name, minimum_useful_peaks = 3, 
    ##     max_iterations = 10, logicle_width = 0.5)
    
    plot(
        exprs(multipeak_results_tf$transformed_data[,"FITC-A"]),
        exprs(multipeak_results_tf$transformed_data[,"Pacific Blue-A"]),
        col=multipeak_results_tf$peak_clusters$cluster, pch='.')

}

\keyword{functions}
