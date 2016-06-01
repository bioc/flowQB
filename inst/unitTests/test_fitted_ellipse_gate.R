test_fitted_ellipse_gate <- function() {
    library('flowCore')
    fcsFilePath <- system.file("extdata", "935289.fcs", package="flowQB")
    myFlowFrame <- read.FCS(fcsFilePath)
    gatedFlowFrame <- fitted_ellipse_gate(myFlowFrame, c('FSC-H', 'SSC-H'))
    checkTrue(description(gatedFlowFrame)$`$TOT` == "7013")
    checkTrue(nrow(gatedFlowFrame) == 7013)
}

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