test_fitted_ellipse_gatek <- function() {
    library('flowCore')
    fcsFilePath <- system.file("extdata", "935289.fcs", package="flowQB")
    myFlowFrame <- read.FCS(fcsFilePath)
    r <- peak_gate(myFlowFrame, 'FSC-H')

    checkTrue(nrow(myFlowFrame[r == TRUE]) == 14133)
    checkTrue(nrow(myFlowFrame[r == FALSE]) == 5867)
}
