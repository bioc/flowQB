# test_fitted_ellipse_gate <- function() {
#     library('flowCore')
#     fcsFilePath <- system.file("extdata", "SSFF_LSRII", "Other_Tests",
#         "933745.fcs", package="flowQBData")
#     myFlowFrame <- read.FCS(fcsFilePath)
#     gatedFlowFrame <- fitted_ellipse_gate(myFlowFrame, c('FSC-H', 'SSC-H'))
#     checkTrue(description(gatedFlowFrame)$`$TOT` == "40526")
#     checkTrue(nrow(gatedFlowFrame) == 40526)
# }

test_fitted_ellipse_gate2 <- function() {
    expected <- "expected"
    tmp <- tryCatch(
        {
            fitted_ellipse_gate("foo", c('FSC-H', 'SSC-H'));
        },
        error = function(ex) {
            expected;
        }
    )
    checkTrue(tmp == expected)
}