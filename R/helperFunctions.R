###############################################################################
## Copyright (c) 2016
## Josef Spidlen, Faysal El Khettabi, Wayne Moore, David Parks, Ryan Brinkman
##
## License
## The software is distributed under the terms of the
## Artistic License 2.0
## http://www.r-project.org/Licenses/Artistic-2.0
## 
## Disclaimer
## This software and documentation come with no warranties of any kind.
## This software is provided "as is" and any express or implied
## warranties, including, but not limited to, the implied warranties of
## merchantability and fitness for a particular purpose are disclaimed.
## In no event shall the copyright holder be liable for any direct,
## indirect, incidental, special, exemplary, or consequential damages
## (including but not limited to, procurement of substitute goods or
## services; loss of use, data or profits; or business interruption)
## however caused and on any theory of liability, whether in contract,
## strict liability, or tort arising in any way out of the use of this
## software.
###############################################################################

## Find density peak in the provided vector of values
find_peak <- function(data, width = 0.5, fraction = 0.1)
{
    ## Additional sanity checks
    if (width >= 1 || width <= 0) 
        stop("The width argument shall be between 0 and 1.")
    if (fraction >= 1 || fraction <= 0)
        stop("The fraction argument shall be between 0 and 1.")
    
    x <- sort(data)
    N <- length(data)
    M <- ceiling(N * fraction)
    M2 <- floor((M + 1) / 2)
    ## Note JS while refactoring: 
    ## Considering that x is sorted, wouldn't first always be 1 after this?
    for (first in 1:(N - M - 1))
        if (x[first] < x[first + 1])
            break
    
    ## Note JS while refactoring: 
    ## Considering that x is sorted, wouldn't last always be N after this?
    for (last in N:(first + M))
        if (x[last - 1] < x[last])
            break
    
    i <- first
    dx <- x[last] - x[first]
    lo <- first
    hi <- last - M
    
    ## Find the densest fraction of cells by minimizing the difference between
    ## the values at the first and "first + fraction size - 1" index of 
    ## sorted data
    for (j in first:(last - M + 1))
    {
        d <- x[j + M - 1] - x[j]
        d
        if (d < dx)
        {
            i <- j
            dx <- d
        }
    }
    
    ## Find the lower bound of the peak based on how much bigger the difference
    ## between the values at the first and "first + fraction size - 1" index of
    ## sorted data is allowed to get. By default, the difference can be
    ## double (i.e., 1 / width)
    for (j in i:first)
    {
        d <- x[j + M - 1] - x[j]
        if (d > dx / width)
        {
            lo <- j
            break
        }
    }
    
    ## Find the upper bound of the peak based on how much bigger the difference
    ## between the values at the first and "first + fraction size - 1" index of
    ## sorted data is allowed to get. By default, the difference can be
    ## double (i.e., 1 / width)
    for (j in i:(last - M + 1))
    {
        d <- x[j + M - 1] - x[j]
        if (d > dx / width)
        {
            hi <- j
            break
        }
    }
    
    list(lo = x[lo + M2], hi = x[hi + M2])
}

## TODO: originally called scatter_gate
## This is applicable to any channels provided as arguments
## to this function, so it has been renamed, but the rest of the code still
## needs to be refactored accordingly
fitted_ellipse_gate <- function(myFlowFrame, channels, R=1)
{
    ## Remove events with negative values in the selected channels 
    ## so that the log transform later doesn't fail
    for (channel in channels)
    {
        values <- exprs(myFlowFrame[,channel])
        myFlowFrame <- myFlowFrame[values>0]
    }
    
    ## Gate by an ellipse (ellipseoid) on the log transformed
    ## values of selected channels
    R2 <- vector(mode="double", nrow(myFlowFrame))
    for (channel in channels)
    {
        values <- log(exprs(myFlowFrame[,channel]))
        fwhm <- find_peak(values)
        center <- (fwhm$hi + fwhm$lo)/2
        radius <- (fwhm$hi - fwhm$lo)/2
        R2 = R2 + ((values - center)/radius)^2
    }

    ## Final gate can be stretched/shrunk by specifying R different from 1
    myFlowFrame <- myFlowFrame[R2 <= R^2]
    ## This does not fix the TOT keyword, so we will take care of that manually
    description(myFlowFrame)$`$TOT` <- as.character(nrow(myFlowFrame))
    myFlowFrame
}
