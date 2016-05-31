is_close_enough <- function(x) {
    if (abs(x[[1]] - x[[2]]) < x[[3]]) return(TRUE)
    FALSE
}

test_calc_mean_sd_duke <- function() {
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
    
    checkTrue(ncol(results) == 36)
    checkTrue(nrow(results) == 4)
    
    expected_rows <- c("total", "scatter gated", "mean", "sd")
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
    
    results['scatter gated',]
    
    expected_scatter_gated <- c(13609, 13612, 13609, 13612, 13610, 2198, 13604, 
        2351, 13605, 2128, 13610, 1700, 13593, 13605, 13609, 1795, 13609, 1937, 
        13608, 2179, 13607, 13608, 13613, 1709, 13606, 13128, 13611, 2295, 
        13609, 5130, 13579, 4945, 13612, 2119, 13606, 13613)
    
    checkTrue(all(results['scatter gated',] == expected_scatter_gated))
    
    expected_mean <- c(
        3.68269806621878, 46.1761248852162, 1.6190191611999,
        32.43507805326, 10.8018460152539, 32.5267045454549, 
        8.19365489811447, 16.5374800637962, 11.157524995952, 
        32.5316901408458, 7.96136292899557, 32.4919117647053, 
        7.92481653176239, 48.9612310519081, 2.22670584941697,
        18.4968684759917, 0.7080594966405, 18.5177304964537, 
        -0.360082650221871, 18.4985673352433, 6.19952875242143,
        33.0154298310049, -0.0542732525403908, 20.4996347699055,
        -0.707918442852837, 67.206587966486, 3.71127478480177,
        16.5666848121936, 1.60514193161841, 18.5416666666671,
        0.954811806334246, 18.5061915592626, 0.854750244019773,
        16.7071302298176, 0.422804537271507, 15.3042879441741
    )

    checkTrue(all(apply(cbind(as.numeric(results['mean',]), 
        expected_mean, 1e-8), 1, is_close_enough)))
    
    expected_sd <- c(
        15.018831382611, 13.3293735563702, 13.3468954509309, 
        12.4763678182567, 36.1962670400232, 0.636028577135425, 
        17.4293656949249, 0.635159693125445, 25.2627769468974, 
        0.635348040182856, 24.6636380580307, 0.638481416076153, 
        27.6728010827464, 13.650047655063, 39.3195268783639, 
        0.637932347466905, 18.6569775140754, 0.644712061644128, 
        23.5381609024289, 0.637546175439074, 10.4918472781137, 
        12.1731829665694, 10.4057245025787, 0.637404309986813, 
        29.5214285585117, 26.4500423101424, 12.5464249628083, 
        0.627783492206857, 10.3801540854037, 2.44696407354404, 
        11.2736179151466, 2.43804911652525, 9.59728636707067, 
        0.896390493635044, 10.7987877822619, 11.2881147959847
    )
    
    checkTrue(all(apply(cbind(as.numeric(results['sd',]), 
        expected_sd, 1e-8), 1, is_close_enough)))

}