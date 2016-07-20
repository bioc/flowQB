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

fit_led <- function(fcs_file_path_list, ignore_channels, 
    dyes, detectors, signal_type, instrument_name, 
    bounds = list(minimum = -100, maximum = 100000),
    minimum_useful_peaks = 3,
    max_iterations = 10, ...)
{
    signal_type <- tolower(signal_type)
    instrument_name <- tolower(instrument_name)
    
    peak.list <- list()
    fluorescences <- c()
    
    for (i in 1:length(fcs_file_path_list))
    {
        fcs <- read.FCS(fcs_file_path_list[[i]])
        peak.list[[i]] <- fcs
    }
    if (length(peak.list) > 0) 
        fluorescences <- pick_parameters(fcs, ignore_channels)
    
    results.list <- get_peak_statistics(peak.list, 
        basename(fcs_file_path_list), fluorescences, bounds, 
        signal_type == "height", instrument_name == "bd accuri", ...)

    bg.result <- data.frame(row.names=c("total", "mean", "sd"))
    for (fluorescence in fluorescences)
    {
        peak.results <- results.list[[fluorescence]]
        bg.result[[fluorescence]] <- c(
            peak.results$N[[1]], peak.results$M[[1]], peak.results$SD[[1]])
    }
    dye.bg.result <- get_results_for_dyes(dyes, detectors, bg.result)

    censored.list <- censor_peak_statistics(
        results.list, fluorescences)
    
    fits <- data.frame(row.names=c(
        "C0", "C0 SE", "C0 P", "C1", "C1 SE", "C1 P", "C2", "C2 SE", "C2 P",
        "C0'", "C0' SE", "C0' P", "C1'", "C1' SE", "C1' P"))
    for (fluorescence in fluorescences)
    {
        results <- results.list[[fluorescence]]
        censored <- censored.list[[fluorescence]]
        Q.R <- vector(mode='double', length=nrow(results))
        L.R <- vector(mode='double', length=nrow(results))
        
        if (usable_rows(censored) >= minimum_useful_peaks)
        {
            ## This is not really needed, but just so that 
            ## R CMD check does not complain about undefined variable
            W <- censored$W

            quadratic.model <- lm(formula = V ~ 1 + M + I(M^2), 
                data = censored, weights = W)
            linear.model <- lm(formula = V ~ 1 + M, 
                data = censored, weights = W)
            q.coef <- coefficients(summary(quadratic.model))
            l.coef <- coefficients(summary(linear.model))
            fits[[fluorescence]] <- c(
                q.coef[1,1], q.coef[1,2], q.coef[1,4], 
                q.coef[2,1], q.coef[2,2], q.coef[2,4], 
                q.coef[3,1], q.coef[3,2], q.coef[3,4],
                l.coef[1,1], l.coef[1,2], l.coef[1,4], 
                l.coef[2,1], l.coef[2,2], l.coef[2,4])

            r <- residuals(quadratic.model)
            for (i in 1:nrow(censored)) 
                if 
                    (censored$Omit[[i]]) Q.R[[i]] <- NA_real_
                else 
                    Q.R[[i]] <- r[[row.names(censored)[[i]]]] * 
                        sqrt(censored$W[[i]])
            
            r <- residuals(linear.model)
            for (i in 1:nrow(censored))
                if (censored$Omit[[i]])
                    L.R[[i]] <- NA_real_
                else
                    L.R[[i]] <- r[[row.names(censored)[[i]]]] * 
                        sqrt(results$W[[i]])
        }
        else
        {
            fits[[fluorescence]][1:nrow(fits)] <- NA_real_
            Q.R[1:length(Q.R)] <- NA_real_
            L.R[1:length(L.R)] <- NA_real_
        }
        
        results$QR <- Q.R
        results$LR <- L.R
        results.list[[fluorescence]] <- results
    }
    
    dye_fits <- get_results_for_dyes(dyes, detectors, fits)

    # Iterated
    iterated_fits <- data.frame(row.names=c(
        "C0", "C0 SE", "C0 P", "C1", "C1 SE", "C1 P", "C2", "C2 SE", "C2 P",
        "C0'", "C0' SE", "C0' P", "C1'", "C1' SE", "C1' P"))
    for (fluorescence in fluorescences)
    {
        results <- results.list[[fluorescence]]
        censored <- censored.list[[fluorescence]]
        Q.R <- vector(mode='double', length=nrow(results))
        L.R <- vector(mode='double', length=nrow(results))
        
        if (usable_rows(censored) >= minimum_useful_peaks)
        {
            ## This is not really needed, but just so that 
            ## R CMD check does not complain about undefined variable
            W <- censored$W

            quadratic.model <- lm(formula = V ~ 1 + M + I(M^2), 
                data = censored, weights = W)
            q.coef <- coefficients(summary(quadratic.model))
            linear.model <- lm(formula = V ~ 1 + M, 
                data = censored, weights = W)
            l.coef <- coefficients(summary(linear.model))

            for (iteration in 1:max_iterations)
            {
                for (i in 1:nrow(censored))
                {
                    V <- q.coef[1,1] + q.coef[2,1] * censored$M[[i]] + 
                        q.coef[3,1] * censored$M[[i]]^2
                    censored$W[[i]] <- (censored$N[[i]] - 1)/(2 * V^2)
                    results$W[[i]] <- censored$W[[i]]
                }
                quadratic.model <- lm(formula = V ~ 1 + M + I(M^2), 
                    data = censored, weights = W)
                qnew.coef <- coefficients(summary(quadratic.model))
                change <- max(abs((qnew.coef[,1] - q.coef[,1])/q.coef[,1]))
                q.coef <- qnew.coef
                if (change < 5E-5) break
            }
            
            for (iteration in 1:max_iterations)
            {
                for (i in 1:nrow(censored))
                {
                    V <- l.coef[1,1] + l.coef[2,1] * censored$M[[i]]
                    censored$W[[i]] <- (censored$N[[i]] - 1)/(2 * V^2)
                    results$W[[i]] <- censored$W[[i]]
                }
                linear.model <- lm(formula = V ~ 1 + M, data = censored, 
                    weights = W)
                lnew.coef <- coefficients(summary(linear.model))
                change <- max(abs((lnew.coef[,1] - l.coef[,1])/l.coef[,1]))
                l.coef <- lnew.coef
                if (change < 5E-5) break
            }

            iterated_fits[[fluorescence]] <- c(
                q.coef[1,1], q.coef[1,2], q.coef[1,4], 
                q.coef[2,1], q.coef[2,2], q.coef[2,4], 
                q.coef[3,1], q.coef[3,2], q.coef[3,4],
                l.coef[1,1], l.coef[1,2], l.coef[1,4], 
                l.coef[2,1], l.coef[2,2], l.coef[2,4])

            r <- residuals(quadratic.model)
            for (i in 1:nrow(censored))
                if (censored$Omit[[i]])
                    Q.R[[i]] <- NA_real_
            else
                Q.R[[i]] <- r[[row.names(censored)[[i]]]] * 
                    sqrt(censored$W[[i]])

            r <- residuals(linear.model)
            for (i in 1:nrow(censored))
                if (censored$Omit[[i]])
                    L.R[[i]] <- NA_real_
            else
                L.R[[i]] <- r[[row.names(censored)[[i]]]] * 
                    sqrt(results$W[[i]])
        }
        else
        {
            iterated_fits[[fluorescence]][1:nrow(iterated_fits)] <- NA_real_
            Q.R[1:length(Q.R)] <- NA_real_
            L.R[1:length(L.R)] <- NA_real_
        }
        
        results['QR-I'] <- Q.R
        results['LR-I'] <- L.R
        results.list[[fluorescence]] <- results
    }
    iterated_dye_fits <- get_results_for_dyes(dyes, detectors, iterated_fits)

    list(
        peak_stats = results.list,
        bg_stats = bg.result,
        dye_bg_stats = dye.bg.result,
        fits = fits,
        dye_fits = dye_fits,
        iterated_fits = iterated_fits,
        iterated_dye_fits = iterated_dye_fits
    )
}

fit_multipeak <- function(fcs_file_path, scatter_channels, ignore_channels,
    N.peaks, dyes, detectors, bounds,
    signal_type, instrument_name,
    minimum_useful_peaks = 3, max_iterations = 10,
    logicle_width = 0.5, ...) {
    if (!file.exists(fcs_file_path)) return()

    signal_type <- tolower(signal_type)
    instrument_name <- tolower(instrument_name)
    
    fcs <- read.FCS(fcs_file_path)
    scatter.gated <- fitted_ellipse_gate(fcs, scatter_channels, 2)
    fluorescences <- pick_parameters(fcs, ignore_channels)

    fluorescence.data <- scatter.gated[,fluorescences]
    logicle = logicleTransform(
        t=parameters(fluorescence.data)$range[[1]],
        w=logicle_width)
    
    x <- exprs(fluorescence.data)
    y <- matrix(logicle(x), nrow=nrow(x), ncol=ncol(x),
        dimnames=list(rownames(x),colnames(x)))
    logicle.data <- new('flowFrame', y, 
        parameters=parameters(fluorescence.data), 
        description=description(fluorescence.data))
    
    peak.list <- list()
    # TODO: 500, 500 could be user-defined parameters?
    ## This seems to result in some warnings that
    ## "Quick-TRANSfer stage steps exceeded maximum (= 3762700)"
    ## but the results seem OK. 
    km <- kmeans(data.frame(y), N.peaks, 500, 500)
    ## I tried rounding up the data a bit as suggested in the kmeans
    ## documentation, but it does not seem to help.
    # km <- kmeans(data.frame(round(y, digits = 3)), N.peaks, 500, 500)
    ## I also tried a different implementation ("Lloyd"), but that comes
    ## with it's own problems and seems to take 5 times as long
    # km <- kmeans(data.frame(y), N.peaks, 500, 500, algorithm="Lloyd")
    for (i in 1:N.peaks)
    {
        peak.list[[i]] <- peak <- fluorescence.data[km$cluster==i]
    }
    
    results.list <- get_peak_statistics(peak.list, 1:N.peaks, fluorescences, 
        bounds, tolower(signal_type) == "height", 
        tolower(instrument_name) == "bd accuri")
    censored.list <- censor_peak_statistics(results.list, fluorescences)

    #fits <- data.frame(row.names=c(
    #    "C0", "C0 SE", "C0 P", "C1", "C1 SE", "C1 P", "C2", "C2 SE", "C2 P",
    #    "C0'", "C0' SE", "C0' P", "C1'", "C1' SE", "C1' P"))
    fits <- data.frame(row.names=c(
        "C0", "C0 SE", "C0 P", "C1", "C1 SE", "C1 P", "C2", "C2 SE", "C2 P"))

    for (fluorescence in fluorescences)
    {
        results <- results.list[[fluorescence]]
        censored <- censored.list[[fluorescence]]
        Q.R <- vector(mode='double', length=nrow(censored))
        # L.R <- vector(mode='double', length=nrow(censored))

        ## This is not really needed, but just so that 
        ## R CMD check does not complain about undefined variable
        W <- censored$W

        if (usable_rows(censored) >= minimum_useful_peaks)
        {
            quadratic.model <- lm(formula = V ~ 1 + M + I(M^2), 
                data = censored, weights = W)
            # linear.model <- lm(formula = V ~ 1 + M, 
            #   data = censored, weights = W)
            q.coef <- coefficients(summary(quadratic.model))
            # l.coef <- coefficients(summary(linear.model))
            
            fits[[fluorescence]] <- c(
                q.coef[1,1], q.coef[1,2], q.coef[1,4], 
                q.coef[2,1], q.coef[2,2], q.coef[2,4], 
                q.coef[3,1], q.coef[3,2], q.coef[3,4]) #,
            #   l.coef[1,1], l.coef[1,2], l.coef[1,4], 
            #   l.coef[2,1], l.coef[2,2], l.coef[2,4])
            
            r <- residuals(quadratic.model)
            for (i in 1:nrow(censored)) 
                if (censored$Omit[[i]]) 
                    Q.R[[i]] <- NA_real_
                else 
                    Q.R[[i]] <- r[[row.names(censored)[[i]]]] * 
                        sqrt(censored$W[[i]])
            
            # r <- residuals(linear.model)
            # for (i in 1:nrow(censored))
            #     if (censored$Omit[[i]]) 
            #         L.R[[i]] <- NA_real_
            #     else
            #         L.R[[i]] <- r[[row.names(censored)[[i]]]] * 
            #             sqrt(results$W[[i]])
            
        }
        else
        {
            fits[[fluorescence]][1:nrow(fits)] <- NA_real_
            Q.R[1:length(Q.R)] <- NA_real_
            # L.R[1:length(L.R)] <- NA_real_
        }
        
        results$QR <- Q.R
        # results$LR <- L.R
        results.list[[fluorescence]] <- results
    }
    
    dye_fits <- get_results_for_dyes(dyes, detectors, fits)

    # Iterated
    # iterated_fits <- data.frame(row.names=c(
    #     "C0", "C0 SE", "C0 P", "C1", "C1 SE", "C1 P", "C2", "C2 SE", "C2 P",
    #     "C0'", "C0' SE", "C0' P", "C1'", "C1' SE", "C1' P"))
    iterated_fits <- data.frame(row.names=c(
        "C0", "C0 SE", "C0 P", "C1", "C1 SE", "C1 P", "C2", "C2 SE", "C2 P"))

    iteration_numbers <- data.frame()
    
    q.iterations = NA
    # l.iterations = NA
    
    for (fluorescence in fluorescences)
    {
        results <- results.list[[fluorescence]]
        censored <- censored.list[[fluorescence]]
        Q.R <- vector(mode='double', length=nrow(censored))
        # L.R <- vector(mode='double', length=nrow(censored))

        if (usable_rows(censored) >= minimum_useful_peaks)
        {
            quadratic.model <- lm(formula = V ~ 1 + M + I(M^2), 
                data = censored, weights = W)
            q.coef <- coefficients(summary(quadratic.model))
            # linear.model <- lm(formula = V ~ 1 + M, 
            #   data = censored, weights = W)
            # l.coef <- coefficients(summary(linear.model))

            for (iteration in 1:max_iterations)
            {
                for (i in 1:nrow(censored))
                {
                    V <- q.coef[1,1] + q.coef[2,1] * censored$M[[i]] + 
                        q.coef[3,1] * censored$M[[i]]^2
                    censored$W[[i]] <- (censored$N[[i]] - 1)/(2 * V^2)
                    results$W[[i]] <- censored$W[[i]]
                }
                quadratic.model <- lm(formula = V ~ 1 + M + I(M^2), 
                    data = censored, weights = W)
                qnew.coef <- coefficients(summary(quadratic.model))
                change <- max(abs((qnew.coef[,1] - q.coef[,1])/q.coef[,1]))
                q.coef <- qnew.coef
                if (change < 5E-5) break
            }
            q.iterations <- iteration
            
            # for (iteration in 1:max_iterations)
            # {
            #     for (i in 1:nrow(censored))
            #     {
            #         V <- l.coef[1,1] + l.coef[2,1] * censored$M[[i]]
            #         censored$W[[i]] <- (censored$N[[i]] - 1)/(2 * V^2)
            #         results$W[[i]] <- censored$W[[i]]
            #     }
            #     linear.model <- lm(formula = V ~ 1 + M, data = censored,
            #         weights = W)
            #     lnew.coef <- coefficients(summary(linear.model))
            #     change <- max(abs((lnew.coef[,1] - l.coef[,1])/l.coef[,1]))
            #     l.coef <- lnew.coef
            #     if (change < 5E-5) break
            # }
            # l.iterations <- iteration

            iterated_fits[[fluorescence]] <- c(
                q.coef[1,1], q.coef[1,2], q.coef[1,4], 
                q.coef[2,1], q.coef[2,2], q.coef[2,4], 
                q.coef[3,1], q.coef[3,2], q.coef[3,4]) #,
            #    l.coef[1,1], l.coef[1,2], l.coef[1,4], 
            #    l.coef[2,1], l.coef[2,2], l.coef[2,4])
            
            r <- residuals(quadratic.model)
            for (i in 1:nrow(censored)) 
                if (censored$Omit[[i]]) 
                    Q.R[[i]] <- NA_real_
                else 
                    Q.R[[i]] <- r[[row.names(censored)[[i]]]] * 
                        sqrt(censored$W[[i]])
            
            # r <- residuals(linear.model)
            # for (i in 1:nrow(censored))
            #     if (censored$Omit[[i]]) 
            #         L.R[[i]] <- NA_real_
            #     else
            #         L.R[[i]] <- r[[row.names(censored)[[i]]]] * 
            #             sqrt(results$W[[i]])
        }
        else
        {
            iterated_fits[[fluorescence]][1:nrow(iterated_fits)] <- NA_real_
            Q.R[1:length(Q.R)] <- NA_real_
            # L.R[1:length(L.R)] <- NA_real_
        }
        
        results['QR-I'] <- Q.R
        # results['LR-I'] <- L.R
        results.list[[fluorescence]] <- results        

        iteration_numbers <- rbind(
            iteration_numbers, 
            data.frame(row.names = c(fluorescence), 
                Q=q.iterations#, L=l.iterations
            ))
    }
    
    iterated_dye_fits <- get_results_for_dyes(dyes, detectors, iterated_fits)
    
    list(
        peak_stats = results.list,
        fits = fits,
        dye_fits = dye_fits,
        iterated_fits = iterated_fits,
        iterated_dye_fits = iterated_dye_fits,
        iteration_numbers = iteration_numbers,
        peak_clusters=km,
        peaks=peak.list, 
        transformed_data=logicle.data
    )
}

fit_spherotech <- function(fcs_file_path, scatter_channels, ignore_channels,
    dyes, detectors, bounds, signal_type, instrument_name, 
    minimum_useful_peaks = 3, max_iterations = 10, logicle_width = 0.5, ...)
{
    fit_multipeak(fcs_file_path, scatter_channels, ignore_channels,
        8, dyes, detectors, bounds, signal_type, instrument_name, 
        minimum_useful_peaks = 3, max_iterations = 10, logicle_width = 0.5, ...)
}

fit_thermo_fischer <- function(fcs_file_path, scatter_channels, ignore_channels,
    dyes, detectors, bounds, signal_type, instrument_name, 
    minimum_useful_peaks = 3, max_iterations = 10, logicle_width = 0.5, ...)
{
    fit_multipeak(fcs_file_path, scatter_channels, ignore_channels,
    6, dyes, detectors, bounds, signal_type, instrument_name, 
    minimum_useful_peaks = 3, max_iterations = 10, logicle_width = 0.5, ...)
}
