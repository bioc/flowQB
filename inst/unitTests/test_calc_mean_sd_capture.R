test_calc_mean_sd_capture_all <- function() {
    library(flowCore)

    file_directory <- system.file("extdata", "example1", 
        "SSFF_LSRII", "SU_2B", package="flowQB")
    fcs_file_path_list <- as.list(file.path(
        file_directory, c("933723.fcs","933725.fcs")))
    scatter_channels_list <- list(c("FSC-A", "SSC-A"), c("FSC-A", "SSC-A"))
    detector_list <- list("APC-A", "APC-Cy7-A")
    dye_list <- list("APC", "APC-Cy7")
    
    results <- calc_mean_sd_capture_all(
        fcs_file_path_list, 
        scatter_channels_list, 
        detector_list, 
        dye_list
    )

    checkTrue(ncol(results) == 2)
    checkTrue(nrow(results) == 8)
    
    expected_rows <- c("total", "scatter gated", "high peak", "low peak",
        "stained mean", "stained sd", "unstained mean", "unstained sd")
    checkTrue(all(expected_rows %in% rownames(results) ==
        rep(TRUE, length(expected_rows))))
    
    expected_cols <- c("APC", "APC-Cy7")
    checkTrue(all(expected_cols %in% colnames(results) ==
        rep(TRUE, length(expected_cols))))
    
    checkTrue(all(results['total',] == 20000))
    
    expected_scatter_gated <- c(15616, 16088)
    checkTrue(all(results['scatter gated',] == expected_scatter_gated))
    
    expected_high_peak <- c(6160, 6837)
    checkTrue(all(results['high peak',] == expected_high_peak))
    
    expected_low_peak <- c(8157, 8253)
    checkTrue(all(results['low peak',] == expected_low_peak))

    expected_stained_mean <- c(41461.08459, 22917.573572)
    checkTrue(all(apply(cbind(as.numeric(results['stained mean',]), 
        expected_stained_mean, 1e-5), 1,
        function(x) {
            if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
        })))
    
    expected_stained_sd <- c(3252.01542, 1438.080887)
    checkTrue(all(apply(cbind(as.numeric(results['stained sd',]), 
        expected_stained_sd, 1e-5), 1,
        function(x) {
            if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
        })))

    expected_unstained_mean <- c(104.26476, 4.109041)
    checkTrue(all(apply(cbind(as.numeric(results['unstained mean',]), 
        expected_unstained_mean, 1e-5), 1,
        function(x) {
            if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
        })))
    
    expected_unstained_sd <- c(23.22988, 33.266498)
    checkTrue(all(apply(cbind(as.numeric(results['unstained sd',]), 
        expected_unstained_sd, 1e-5), 1,
        function(x) {
            if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
        })))
}

test_calc_mean_sd_capture <- function() {
    library(flowCore)

    fcs_file_path <- system.file("extdata", "example1", "SSFF_LSRII", "SU_2B",
        "933723.fcs", package="flowQB")

    scatter_channels <- c("FSC-A", "SSC-A")
    detector <- "APC-A"
    dye <- "APC"
    
    results <- calc_mean_sd_capture(
        fcs_file_path, scatter_channels, detector, dye)
    
    checkTrue(colnames(results) == "APC")
    checkTrue(length(colnames(results)) == 1)
    checkTrue(results['total',] == 20000)
    checkTrue(results['high peak',] == 6160)
    checkTrue(results['low peak',] == 8157)
    
    checkTrue(results['stained mean',] >= 41461.0845 && 
        results['stained mean',] <= 41461.0846)
    checkTrue(results['stained sd',] >= 3252.0154 && 
        results['stained sd',] <= 3252.0155)
    checkTrue(results['unstained mean',] >= 104.2647 && 
        results['unstained mean',] <= 104.2648)
    checkTrue(results['unstained sd',] >= 23.2298 && 
        results['unstained sd',] <= 23.2299)
}