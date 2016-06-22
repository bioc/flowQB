test_calc_mean_sd_197 <- function() {
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
    
    test_column <- 14
    test_row <- 14
    folder <- xlsx[[test_column]][[test_row]]
    beads_file_name <- xlsx[[test_column]][[test_row+1]]
    scatter_channels <- c(
        xlsx[[test_column]][[test_row+2]], 
        xlsx[[test_column]][[test_row+3]])
    
    fcs_path <- system.file("extdata", instrument_folder, folder,
        beads_file_name, package="flowQBData")
    
    results <- calc_mean_sd_duke(fcs_path, scatter_channels, ignore_channels)
    
    checkTrue(ncol(results) == 36)
    checkTrue(nrow(results) == 5)
    
    expected_rows <- c("total", "scatter gated", "peak gated", "mean", "sd")
    checkTrue(all(expected_rows %in% rownames(results) == 
        rep(TRUE, length(expected_rows))))
    
    expected_cols <- c("FITC-A", "FITC-H", "PerCP-Cy5-5-A", 
        "PerCP-Cy5-5-H", "Pacific Blue-A", "Pacific Blue-H", 
        "Aqua Amine-A", "Aqua Amine-H", "Pacific Orange-A", 
        "Pacific Orange-H", "QDot 585-A", "QDot 585-H", "QDot 605-A", 
        "QDot 605-H", "QDot 655-A", "QDot 655-H", "QDot 705-A", 
        "QDot 705-H", "QDot 800-A", "QDot 800-H", "APC-A", "APC-H", 
        "APC-Cy5-5-A", "APC-Cy5-5-H", "APC-Cy7-A", "APC-Cy7-H", "PE-A", 
        "PE-H", "PE-Texas-Red-A", "PE-Texas-Red-H", "PE-Cy5-A", 
        "PE-Cy5-H", "PE-Cy5-5-A", "PE-Cy5-5-H", "PE-Cy7-A", "PE-Cy7-H")
    
    checkTrue(all(expected_cols %in% colnames(results) == 
        rep(TRUE, length(expected_cols))))

    checkTrue(all(results['total',] == 20000))
    checkTrue(all(results['scatter gated',] == 14368))

    expected_peak_gated <- c(
        14350, 14360, 14367, 14368, 14300, 14313, 14341, 14352, 14355, 14362, 
        14359, 14360, 14366, 14366, 14368, 14368, 14368, 14368, 14367, 14367, 
        14236, 14281, 14223, 14246, 14284, 14326, 14364, 14367, 14366, 14367, 
        14362, 14362, 14363, 14366, 14368, 14368)
    checkTrue(all(results['peak gated',] == expected_peak_gated))

    expected_mean <- c(
        6370.36831790505, 6377.23903203361, 7101.31971088836,
        7277.01487473904, 94005.67182105, 93823.4342852161,
        34835.699258245, 35052.4596760147, 28061.0180951038,
        27183.213054831, 19621.749678871, 19162.6252611421,
        13021.7946194924, 12764.8933356534, 3892.29947296548,
        3815.98329853862, 2118.01693350752, 2091.45746346555,
        732.409636975949, 731.139103958238, 9982.8642829339,
        9664.0848140044, 4530.85891087749, 4403.09914020006,
        3507.15254813337, 3520.85447565871, 61988.38277443,
        62454.7738147017, 28851.4338472791, 29303.0372335797,
        16070.1756378679, 16433.1580504778, 10269.9951403034,
        10469.3153819384, 3389.87192855549, 3480.45807237303
    )

    checkTrue(all(apply(cbind(as.numeric(results['mean',]), 
        expected_mean, 1e-6), 1, function(x) {
            if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
        })))
    
    expected_sd <- c(
        111.759086731826, 140.080599596617, 211.22202829812,
        279.377205174258, 699.615686700079, 773.950353049118,
        381.395955069773, 418.89862222668, 380.20348969832,
        396.138661334339, 314.231399263485, 327.618956840525,
        277.841065786121, 288.720402622319, 173.154703595224,
        174.521682841442, 116.074183373492, 120.154461299576,
        105.303421569486, 106.398899915912, 593.787779173592,
        638.275704605673, 282.957314812784, 302.502816606589,
        298.159991218402, 363.888951175724, 1047.89984052113,
        1130.00478776555, 489.277732577995, 531.32377845144,
        334.649736189545, 390.483492959047, 211.627813695109,
        247.740985058674, 117.529599386806, 154.644307939885
    )
    
    checkTrue(all(apply(cbind(as.numeric(results['sd',]), 
        expected_sd, 1e-6), 1, function(x) {
            if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
        })))

}