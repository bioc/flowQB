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

get_results_for_dyes <- function(dyes, detectors, results)
{
    dye.results <- data.frame(row.names=row.names(results))
    for (i in 1:min(length(dyes),length(detectors))) 
    {
        dye <- dyes[i]
        detector <- detectors[i]
        if (is.na(detector)) {
            dye.results[[dye]] <- vector(mode='numeric', nrow(results))
        } else {
            dye.results[[dye]] <- results[[detector]]
        }
    }
    dye.results
}

# TODO 
# peak.list => flowFrame_list
# peak.names => file_names
# maximum.cv.area => maximum_cv_area
# maximum.cv.height => maximum_cv_height
get_peak_statistics <- function(peak.list, peak.names = NULL, parameters, 
    bounds, is.height = FALSE, is.accuri = FALSE,
    maximum.cv.area = .65, maximum.cv.height = .65)
{
    if (is.null(peak.names)) {
        for (i in 1:length(peak.list)) 
            peak.names <- c(peak.names, peak.list[[i]]@description['$FIL'])
    }

    results.list <- list()
    for (fluorescence in parameters)
        results.list[[fluorescence]] <- data.frame(row.names=peak.names)
    
    for (i in 1:length(peak.list))
    {
        peak <- peak.list[[i]]
        for (fluorescence in parameters)
        {
            data <- exprs(peak[,fluorescence])
            results.list[[fluorescence]]$N[[i]] <- N <- nrow(data)
            out <- getOutliers(as.vector(data), distribution="normal")
            results.list[[fluorescence]]$M[[i]] <- M <- out$mu
            results.list[[fluorescence]]$SD[[i]] <- SD <- out$sigma
            results.list[[fluorescence]]$V[[i]] <- V <- SD^2
            results.list[[fluorescence]]$W[[i]] <- (N - 1)/(2 * V^2)
            results.list[[fluorescence]]$Omit[[i]] <- FALSE
        }
    }
    
    for (fluorescence in parameters)
    {
        results <- results.list[[fluorescence]]
        results.list[[fluorescence]] <- results[order(results$M),]
    }
    
    for (fluorescence in parameters)
    {
        results <- results.list[[fluorescence]]
        ## TODO
        ## Seems like the next 2 if statements could be combined to 
        ## something like
        ## if (is.height || length(grep("-H", fluorescence, fixed=TRUE)) > 0 ||
        ## is.accuri)
        if ((!is.height && length(grep("-H", fluorescence, fixed=TRUE)) > 0) 
            || is.accuri)
        {
            lowest.mean <- results$M[[1]]
            for (i in 1:nrow(results))
            {
                this.mean <- results$M[[i]]
                if (this.mean < 10 * lowest.mean) results$Omit[[i]] <- TRUE
            }
        }
        if (is.height)
        {
            lowest.mean <- results$M[[1]]
            for (i in 1:nrow(results))
            {
                this.mean <- results$M[[i]]
                if (this.mean < 10 * lowest.mean) results$Omit[[i]] <- TRUE
            }
        }
        for (i in 1:nrow(results))
        {
            if (is.height) 
                maximum.cv <- maximum.cv.height 
            else 
                maximum.cv <- maximum.cv.area
            omit <- results$Omit[[i]]
            if (!omit && (results$M[[i]] > bounds$maximum || 
                results$SD[[i]] == 0)) omit <- TRUE
            if (!omit && results$M[[i]] < bounds$minimum) 
                omit <- TRUE else CV <- results$SD[[i]]/results$M[[i]]
            if (!omit && (is.height || is.accuri || 
                length(grep("-H", fluorescence, fixed=TRUE)) > 0) && 
                CV > maximum.cv) omit <- TRUE
            if (omit) results$Omit[[i]] <- TRUE
        }
        results.list[[fluorescence]] <- results
    }
    
    results.list
}

censor_peak_statistics <- function(results.list, parameters)
{
    censored.list <- list()
    for (fluorescence in parameters)
    {
        results <- results.list[[fluorescence]]
        censored <- data.frame(results, check.rows=FALSE, check.names=FALSE)
        for (i in 1:nrow(censored))
        {
            if(censored$Omit[[i]]) censored$M[[i]] <- censored$SD[[i]] <- 
                censored$V[[i]] <- censored$W[[i]] <- NA_real_
        }
        censored.list[[fluorescence]] <- censored
    }
    censored.list
}

usable_rows <- function(results)
{
    count <- 0
    for (i in 1:nrow(results)) if (!is.na(results$M[[i]])) count <- count + 1
    count
}

qb_from_fits <- function(fits) {
    if (length(rownames(fits)) >= 7 && 
        (all(rownames(fits)[c(1,4,7)] == c("c0", "c1", "c2")))) {
        sapply(fits, function(x) { 
            if(length(x) >= 15) {
                list(
                    'q_QI'    = 1 / x[4],        # QI    = 1/c1
                    'q_BSpe'  = x[1] / (x[4]^2), # BSpe  = c0/c1^2
                    'q_CV0sq' = x[7],            # CV0^2 = c2
                    'l_QI'    = 1 / x[13],
                    'l_BSpe'  = x[10] / (x[13]^2)
                )
            } else {
                list(
                    'q_QI'    = 1 / x[4],        # QI    = 1/c1
                    'q_BSpe'  = x[1] / (x[4]^2), # BSpe  = c0/c1^2
                    'q_CV0sq' = x[7]             # CV0^2 = c2
                )
            }
        })
    } else {
        stop(paste("Incompatible fits (the provided fits are not compatible",
            "with the return values of the fit_led or fit_multipeak functions"))
    }
}
