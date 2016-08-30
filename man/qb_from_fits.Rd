\name{qb_from_fits}
\alias{qb_from_fits}

\title{ 
    Extract Q, B and the intrinsic CV0 from fitting results
}

\description{
    This function can be used to compute 
    (1) flow cytometer's detection efficiency (i.e., Q, the statistical number 
    of single photoelectrons (Spe) generated per unit of dye in the sample),
    (2) background illumination (i.e., B, the background light in dye 
    equivalents that sets the minimum variance that underlies all 
    measurements), and 
    (3) the intrinsic CV0 (i.e, variance in the signal produced by the
    variation in dye amount of beads containing a "fixed" level of dye plus the
    illumination variations due to particles taking different flow paths
    through the laser beam)
    from fitting results produced by either LED data fitting (\code{fit_led} 
    function) or bead data fitting (\code{fit_beads}, \code{fit_spherotech}
    or \code{fit_thermo_fisher} functions). One can calculate based on any
    of the fits (i.e, the \code{fits}, \code{dye_fits}, \code{iterated_fits}
    or \code{iterated_dye_fits} items from the result list).
}

\usage{
    qb_from_fits(fits)
}

\arguments{
    \item{fits}{
Fitting results as produced by either LED data fitting (\code{fit_led} 
function) or bead data fitting (\code{fit_beads}, \code{fit_spherotech}
or \code{fit_thermo_fisher} functions). One can calculate based on any
of the fits provided as part of the results of these functions, namely
the \code{fits}, \code{dye_fits}, \code{iterated_fits} or 
\code{iterated_dye_fits} items from the result list.
    }
}

\details{
    As explained in our paper, 
    this method calculats QI as 1/c1, BSpe  = c0/c1^2 and CV0^2 = c2
    from the results of quadratic fit for the measured means and variances 
    to a statistical model involving the Poisson distribution expectations 
    for the relation between them.
    For fitting results containing both, a quadratic and a linear fit
    coeficients, (i.e, results of LED fitting), this method also includes
    QI and BSpe from the linear model.
}

\value{
    The value is a matrix with columns corresponding to the columns of the
    input data frame (i.e., the names of the dyes). The rows are as follows:
    q_QI as the QI from the fits of the quadratic model,
    q_BSpe as the BSpe from the fits of the quadratic model,
    q_CV0sq as the CV0^2 from the fits of the quadratic model,
    l_QI as the QI from the fits of the linear model (provided only if linear
    fit coefficients are also present in the input data frame, i.e., for LED
    fits only), and
    l_BSpe as the BSpe from the fits of the linear model (again, provided only
    if linear fit coefficients are also present in the input data frame, 
    i.e., for LED fits only)
}

\author{
    Josef Spidlen, Wayne Moore, Faysal El Khettabi 
}

\seealso{
    \code{\link{fit_led}},
    \code{\link{fit_beads}},
    \code{\link{fit_spherotech}},
    \code{\link{fit_thermo_fisher}}
}

\examples{
## Example is based on LED data from the flowQBData package
fcs_directory <- system.file("extdata", "SSFF_LSRII", "LED_Series", 
    package="flowQBData")
fcs_file_path_list <- list.files(fcs_directory, "*.fcs", full.names= TRUE)
## We are working with these FCS files:
basename(fcs_file_path_list)

## Various house keeping information
## - Which channels should be ignored, typically the non-fluorescence
##    channels, such as the time and the scatter channels
ignore_channels <- c("Time", 
    "FSC-A", "FSC-W", "FSC-H", 
    "SSC-A", "SSC-W", "SSC-H")
## - Which dyes would you typically use with the detectors
dyes <- c("APC", "APC-Cy7", "APC-H7", "FITC", "PE", "PE-Cy7", "PerCP", 
    "PerCP-Cy55", "V450", "V500-C")
## - What are the corresponding detectors, provide a vector of short channel 
## names, i.e., values of the $PnN FCS keywords.
detectors <- c("APC-A", "APC-Cy7-A", "APC-Cy7-A", "FITC-A", "PE-A", "PE-Cy7-A",
    "PerCP-Cy5-5-A", "PerCP-Cy5-5-A", "Pacific Blue-A", "Aqua Amine-A")
## - The signal type that you are looking at (Area, Height or Width)
signal_type <- "Area"
## - The instrument make/model
instrument_name <- 'LSRII'
## - Set the minimum and maximum values, peaks with mean outsize of this range
## will be ignored
bounds <- list(minimum = -100, maximum = 100000)
## - The minimum number of usable peaks (represented by different FCS files
## in case of an LED pulser) required in order for a fluorescence channel 
## to be included in the fitting. Peaks with mean expression outside of the 
## bounds specified above are omitted and therefore not considered useful
minimum_fcs_files <- 3 # The default 3 seems to be work well in typical cases
## - What is the maximum number of iterations for iterative fitting with
## weight adjustments
max_iterations <- 10 # The default 10 seems to be enough in typical cases

## Now, let's calculate the fitting
led_results <- fit_led(fcs_file_path_list, ignore_channels, dyes,
    detectors, signal_type, instrument_name, bounds = bounds,
    minimum_useful_peaks = minimum_fcs_files, max_iterations = max_iterations)

qb_from_fits(led_results$iterated_dye_fits)
}

\keyword{functions}
