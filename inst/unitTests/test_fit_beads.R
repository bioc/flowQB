# test_fit_beads <- function() {
#     library(flowCore)
#     library(flowQBData)
#     
#     fcs_file_path <- system.file("extdata", "SSFF_LSRII", "Other_Tests", 
#         "933745.fcs", package="flowQBData")
#     scatter_channels <- c("FSC-A", "SSC-A")
#     ignore_channels <- c("Time", "FSC-A", "FSC-W", "FSC-H", 
#         "SSC-A", "SSC-W", "SSC-H")
#     dyes <- c("APC", "APC-Cy7", "APC-H7", "FITC", "PE", "PE-Cy7", "PerCP",
#         "PerCP-Cy55", "V450", "V500-C")
#     detectors <- c("APC-A", "APC-Cy7-A", "APC-Cy7-A", "FITC-A", "PE-A", 
#         "PE-Cy7-A", "PerCP-Cy5-5-A", "PerCP-Cy5-5-A", "Pacific Blue-A",
#         "Aqua Amine-A")
#     bounds <- list(minimum = -100, maximum = 100000)
#     signal_type <- "Area"
#     instrument_name <- 'LSRII'
#     
#     set.seed(123456)
#     beads_results <- fit_spherotech(fcs_file_path, scatter_channels, 
#         ignore_channels, dyes, detectors, bounds, 
#         signal_type, instrument_name, minimum_useful_peaks = 3, 
#         max_iterations = 10, logicle_width = 0.5)
#     
#     ## This would be the same thing except it gives the option of specifying
#     ## the number of peaks
#     ## beads_results <- fit_beads(fcs_file_path, scatter_channels, 
#     ##     ignore_channels, 8, dyes, detectors, bounds, 
#     ##     signal_type, instrument_name, minimum_useful_peaks = 3, 
#     ##     max_iterations = 10, logicle_width = 0.5)
# 
#     checkTrue(
#         apply(cbind(
#             sum(beads_results$fits, na.rm=TRUE), 9352.36544890927, 1e-5), 
#             1, function(x) {
#                 if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
#             }
#         ))
#     
#     checkTrue(
#         apply(cbind(
#             sum(beads_results$dye_fits, na.rm=TRUE), 
#             4849.85545771564, 1e-5),
#             1, function(x) {
#                 if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
#             }
#         ))
#     
#     checkTrue(
#         apply(cbind(
#             sum(beads_results$iterated_fits, na.rm=TRUE), 
#             9227.50586712488, 1e-5),
#             1, function(x) {
#                 if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
#             }
#         ))
#     
#     checkTrue(
#         apply(cbind(
#             sum(beads_results$iterated_dye_fits, na.rm=TRUE), 
#             4886.04514679602, 1e-5), 
#             1, function(x) {
#                 if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
#             }
#         ))
#     
#     checkTrue(
#         apply(cbind(
#             sum(beads_results$iteration_numbers, na.rm=TRUE), 
#             119, 1e-5), 
#             1, function(x) {
#                 if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
#             }
#         ))
#     
# 
#     checkTrue(
#         apply(cbind(
#             sum(beads_results$peak_clusters$cluster, na.rm=TRUE), 
#             338927, 1e-5), 
#             1, function(x) {
#                 if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
#             }
#         ))
#     
# 
#     checkTrue(
#         apply(cbind(
#             length(beads_results$peaks), 
#             8, 1e-5), 
#             1, function(x) {
#                 if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
#             }
#         ))
#     
#     
#     checkTrue(
#         apply(cbind(
#             nrow(beads_results$peaks[[1]]), 
#             9929, 1e-5), 
#             1, function(x) {
#                 if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
#             }
#         ))
#     
#     checkTrue(
#         apply(cbind(
#             sum(beads_results$transformed_data@exprs), 
#             6072572.6847956, 1e-5), 
#             1, function(x) {
#                 if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
#             }
#         ))
#     
#     ## Thermo-Fisher test:
#     fcs_file_path <- system.file("extdata", "SSFF_LSRII", "Other_Tests", 
#         "933747.fcs", package="flowQBData")
#     
#     beads_results_tf <- fit_thermo_fisher(fcs_file_path, scatter_channels, 
#         ignore_channels, dyes, detectors, bounds, 
#         signal_type, instrument_name, minimum_useful_peaks = 3,
#         max_iterations = 10, logicle_width = 0.5)
#     
#     ## This would be the same thing except it gives the option of specifying
#     ## the number of peaks
#     ## beads_results_tf <- fit_beads(fcs_file_path, scatter_channels, 
#     ##     ignore_channels, 6, dyes, detectors, bounds, 
#     ##     signal_type, instrument_name, minimum_useful_peaks = 3,
#     ##     max_iterations = 10, logicle_width = 0.5)
#  
#     checkTrue(
#         apply(cbind(
#             sum(beads_results_tf$fits, na.rm=TRUE), 16706.6395944025, 1e-5),
#             1, function(x) {
#                 if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
#             }
#         ))
#     
#     checkTrue(
#         apply(cbind(
#             sum(unlist(qb_from_fits(beads_results$iterated_dye_fits))),
#             624.505385861861,
#             1e-5), 1, function(x) {
#                 if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
#             }
#         ))
#     
#     checkTrue(
#         apply(cbind(
#             sum(unlist(qb_from_fits(beads_results$dye_fits))),
#             611.787410298806,
#             1e-5), 1, function(x) {
#                 if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
#             }
#         ))
# 
#     checkTrue(
#         apply(cbind(
#             sum(unlist(qb_from_fits(beads_results$iterated_fits))),
#             1064.56472334922,
#             1e-5), 1, function(x) {
#                 if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
#             }
#         ))
# 
#     checkTrue(
#         apply(cbind(
#             sum(unlist(qb_from_fits(beads_results$fits))),
#             1058.67952214706,
#             1e-5), 1, function(x) {
#                 if (abs(x[[1]] - x[[2]]) < x[[3]]) TRUE else FALSE
#             }
#         ))
#     
#     ## TODO 
#     ## Consider adding more checks for the Thermo-Fisher results?
# 
# }
# 
# test_fit_beads_na <- function() {
#     checkTrue(is.null(
#         fit_beads(
#             fcs_file_path = "nonExistantFile.xxx", 
#             scatter_channels = NA, 
#             ignore_channels = NA,
#             N_peaks = 10,
#             dyes = NA,
#             detectors = NA, 
#             bounds = NA,
#             signal_type = "area", 
#             instrument_name = "aria",
#             minimum_useful_peaks = 3, 
#             max_iterations = 10,
#             logicle_width = 1.0
#         )
#     ))
# }