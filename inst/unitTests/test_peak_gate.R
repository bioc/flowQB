# test_peak_gate <- function() {
#     library('flowCore')
#     fcsFilePath <- system.file("extdata", "SSFF_LSRII", "Other_Tests",
#         "933745.fcs", package="flowQBData")
#     myFlowFrame <- read.FCS(fcsFilePath)
# 
#     r <- peak_gate(myFlowFrame, 'FSC-H')
#     checkTrue(nrow(myFlowFrame[r == TRUE]) == 66144)
#     checkTrue(nrow(myFlowFrame[r == FALSE]) == 33856)
# 
#     r2 <- peak_gate(myFlowFrame, 'FSC-H', 2)
#     checkTrue(nrow(myFlowFrame[r2 == TRUE]) == 86810)
#     checkTrue(nrow(myFlowFrame[r2 == FALSE]) == 13190)
# }
# 
# test_peak_gate2 <- function() {
#     expected <- "expected"
#     tmp <- tryCatch(
#         {
#             peak_gate("foo", 'FSC-H');
#         },
#         error = function(ex) {
#             expected;
#         }
#     )
#     checkTrue(tmp == expected)
# }