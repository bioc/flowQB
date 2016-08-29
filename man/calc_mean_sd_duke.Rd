\name{calc_mean_sd_duke}
\alias{calc_mean_sd_duke}

\title{ 
    Calculate the mean and the standard deviation of calibration beads 
    from Duke.
}

\description{
    Currently, this is the same calculation as \code{calc_mean_sd_197}.
}

\usage{
    calc_mean_sd_duke(fcs_file_path, scatter_channels, ignore_channels)
}

\arguments{
    \item{fcs_file_path}{
A character string specifying the file path
to the FCS file with the acquired bead data.    
    }
    \item{scatter_channels}{
A vector of 2 short channel names (values of
the $PnN keywords) specifying the 2 channels that should not be used to gate
the main bead population. The first channel should be a forward scatter
channel, the second one should be a side scatter channel.
    }
    \item{ignore_channels}{
A vector of short channel names (values of
the $PnN keywords) specifying channels that should not be considered for
the fitting procedure. Normally, those should be all non-fluorescence channels,
such as the time and the (forward and side) scatter channels.
    }
}

\details{
    Currently, this is the same calculation as \code{calc_mean_sd_197}.
}

\value{
    The result is a data frame with columns corresponding to short channel
    names of channels from the input FCS file except those specified by the
    ignore_channels parameter. The rows include the total number of events,
    the number of events in the FSC/SSC ellipse gate, the number of 
    events in the peak gate (which can vary slighly among the different 
    channels), the mean and the standard deviation.
}

\author{
    Wayne Moore, Faysal El Khettabi, Josef Spidlen
}

\seealso{
    \code{\link{fitted_ellipse_gate}},
    \code{\link{calc_mean_sd_197}}
}

\examples{
    library(flowCore)
    library(xlsx)
    ## library(flowQBData)
    ## flowQBData is available since BioConductor 3.4, please install it
    ## manually in order to be able to these examples

#   inst_xlsx_path <- system.file("extdata", 
#       "140126_InstEval_Stanford_LSRIIA2.xlsx", package="flowQBData")
#   xlsx <- read.xlsx(inst_xlsx_path, 1, headers=FALSE, stringsAsFactors=FALSE)
#    
#   ignore_channels_row <- 9
#   ignore_channels <- vector()
#   i <- 1
#   while(!is.na(xlsx[[i+4]][[ignore_channels_row]])) {
#       ignore_channels[[i]] <- xlsx[[i+4]][[ignore_channels_row]]
#       i <- i + 1
#   }
#    
#   instrument_folder_row <- 9
#   instrument_folder_col <- 2
#   instrument_folder <- xlsx[[instrument_folder_col]][[instrument_folder_row]]
#
#   test_column <- 13
#   test_row <- 14
#   folder <- xlsx[[test_column]][[test_row]]
#   beads_file_name <- xlsx[[test_column]][[test_row+1]]
#   scatter_channels <- c(
#       xlsx[[test_column]][[test_row+2]], 
#       xlsx[[test_column]][[test_row+3]])
#
#   fcs_path <- system.file("extdata", instrument_folder, folder,
#       beads_file_name, package="flowQBData")
#
#   results <- calc_mean_sd_duke(fcs_path, scatter_channels, ignore_channels)
}

\keyword{functions}
