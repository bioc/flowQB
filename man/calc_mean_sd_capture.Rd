\name{calc_mean_sd_capture}
\alias{calc_mean_sd_capture}

\title{ 
    Calculate the mean and the standard deviation 
    for the stained and unstained population
    of a specified channel in an FCS file.
}

\description{
    This function calculates the mean and the standard deviation of
    two populations - one unstaned (low) and one stained (high) for
    a specified channel (specified by the value of the detector parameter).
    These populations are derived by first fitting an ellipse gate on the
    2 specified scatter channels and then splitting for low and high
    based on the specified detector.
}

\usage{
    calc_mean_sd_capture(fcs_file_path, scatter_channels, detector, dye)
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
    \item{detector}{
A character string specifying which channel to split for the low/high gate, 
which is also the channel that we calculate the mean and the standard deviation
for.
    }
    \item{dye}{
A character string specifying the desired column heading of the result.
    }
}

\details{
    This function first fits an ellipse gate on the 2 specified scatter
    channels. This scatter gated population is then split to 2 (low=unstained and
    high=stained) in each of the fluorescence channels, and a peak 
    gate is applied in order to isolate the high and the low peaks.
    Finally, the
    getOutliers method from the extremevalues package is used in order to
    calculate the mean and the standard deviation of both the
    stained and unstained population of each of the fluorescence channels.
    The value of the detector argument is used to determine which channel
    to work with, the value of the dye argument is only used to specify
    the column name of the result. These shall correspond to each other, e.g.,
    "FITC-A" may be the value of the detector, which shall correspond to the
    short channel name of a channel in the FCS file. FITC may then be the
    desidered dye name, which will end be used for the column heading.
}

\value{
    The result is a data frame with a single column, the heading of the
    column corresponds to the value of the dye argument. The rows include the
    total number of events,  the number of events in the FSC/SSC ellipse gate,
    the number of events in the high peak gate and low peak gate, the stained
    mean and stained standard deviation (based on the high peak gate), and
    finally the unstained mean and unstained standard deviation (based on
    the low peak gate).
}

\author{
    Wayne Moore, Faysal El Khettabi, Josef Spidlen
}

\seealso{
    \code{\link{calc_mean_sd_capture_all}}
}

\examples{
    library(flowCore)
    library(flowQBData)

    fcs_file_path <- system.file("extdata", "SSFF_LSRII", "SU_2B",
        "933723.fcs", package="flowQBData")

    scatter_channels <- c("FSC-A", "SSC-A")
    detector <- "APC-A"
    dye <- "APC"
    
    results <- calc_mean_sd_capture(
        fcs_file_path, scatter_channels, detector, dye)
}

\keyword{functions}

