calc_mean_sd_background <- function() {
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
    
    test_column <- 15
    test_row <- 14
    folder <- xlsx[[test_column]][[test_row]]
    file_name <- xlsx[[test_column]][[test_row+1]]

    fcs_path <- system.file("extdata", "example1", 
        instrument_folder, folder, file_name, package="flowQB")
    
    results <- calc_mean_sd_background(fcs_path, ignore_channels)
    
    checkTrue(ncol(results) == 36)
    checkTrue(nrow(results) == 3)
    
    expected_rows <- c("total", "mean", "sd")
    checkTrue(all(expected_rows %in% rownames(results) == 
        rep(TRUE, length(expected_rows))))
    
    expected_cols <- c(
        "FITC-A", "FITC-H", "PerCP-Cy5-5-A", "PerCP-Cy5-5-H", "Pacific Blue-A",
        "Pacific Blue-H", "Aqua Amine-A", "Aqua Amine-H", "Pacific Orange-A",
        "Pacific Orange-H", "QDot 585-A", "QDot 585-H", "QDot 605-A",
        "QDot 605-H", "QDot 655-A", "QDot 655-H", "QDot 705-A", "QDot 705-H",
        "QDot 800-A", "QDot 800-H", "APC-A", "APC-H", "APC-Cy5-5-A",
        "APC-Cy5-5-H", "APC-Cy7-A", "APC-Cy7-H", "PE-A", "PE-H",
        "PE-Texas-Red-A", "PE-Texas-Red-H", "PE-Cy5-A", "PE-Cy5-H",
        "PE-Cy5-5-A", "PE-Cy5-5-H", "PE-Cy7-A", "PE-Cy7-H"
    )
    
    checkTrue(all(expected_cols %in% colnames(results) == 
        rep(TRUE, length(expected_cols))))

    checkTrue(all(results['total',] == 20000))

    expected_mean <- c(
        -0.774549990139915, 39.3061249999997, -1.29128123424947,
        28.7233749999996, -0.64808248580985, 23.644624999999,
        -1.5814800118731, 20.7499375000008, -1.49790373907237,
        35.0721874999956, -0.730796253167143, 34.4681250000003,
        -1.01593874091661, 44.9691249999991, -2.26630497633682,
        34.9159374999996, -1.30530749563886, 22.3330000000003,
        -1.8858562432156, 27.8691874999991, -0.464514368455832,
        18.6968749999985, -0.359954994838687, 22.9634999999997,
        -1.27972123260054, 68.7965625000004, -0.810117518872044,
        29.8804374999999, -0.254541887432349, 20.6448125000006,
        -0.461268759161192, 24.0135000000003, -0.75586940619344,
        21.0438750000007, -0.603208141416327, 23.7435000000002
    )

    checkTrue(all(apply(cbind(as.numeric(results['mean',]), 
        expected_mean, 1e-5), 1, function(x) {
            if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
        })))
    
    expected_sd <- c(
        10.4514632846542, 11.8778481273538, 11.1050418892498,
        10.540268311226, 33.6198862243458, 8.66119443397132,
        15.8149813011837, 8.49131773022434, 20.6605489391602,
        10.4387045775763, 19.2064528216041, 11.9113692921393,
        22.3738018182496, 14.2354888254652, 37.2405813554801,
        12.5898174791433, 17.4662741073315, 6.81466799907188,
        22.3803127207349, 10.6926381091479, 9.30397135933276,
        12.1099794185227, 9.02058351108572, 6.82537682979835,
        25.0498609743255, 27.6171853891848, 10.134432704732,
        10.9571135623987, 8.7266230817367, 8.28496617042062,
        9.87500041491691, 6.54897886664431, 8.51897027897363,
        7.95449635551379, 10.2552401284315, 6.12760137403496
    )
    
    checkTrue(all(apply(cbind(as.numeric(results['sd',]), 
        expected_sd, 1e-5), 1, function(x) {
            if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
        })))

}