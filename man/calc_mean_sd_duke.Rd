\name{calc_mean_sd_duke}
\alias{calc_mean_sd_duke}

\title{ 
    TODO
}

\description{
    TODO.
}

\usage{
    calc_mean_sd_duke(fcs_file_path, scatter_channels, ignore_channels)
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
    \code{\link{fitted_ellipse_gate}}
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

    test_column <- 13
    test_row <- 14
    folder <- xlsx[[test_column]][[test_row]]
    beads_file_name <- xlsx[[test_column]][[test_row+1]]
    scatter_channels <- c(
        xlsx[[test_column]][[test_row+2]], 
        xlsx[[test_column]][[test_row+3]])

    fcs_path <- system.file("extdata", "example1", 
        instrument_folder, folder, beads_file_name, package="flowQB")

    results <- calc_mean_sd_duke(fcs_path, scatter_channels, ignore_channels)
}

\keyword{functions}
