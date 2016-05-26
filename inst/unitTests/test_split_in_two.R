test_fitted_ellipse_gatek <- function() {
    library('flowCore')
    fcsFilePath <- system.file("extdata", "935289.fcs", package="flowQB")
    myFlowFrame <- read.FCS(fcsFilePath)
    r1 <- split_in_two(myFlowFrame, 'FSC-H')
    r2 <- split_in_two(exprs(myFlowFrame[,'SSC-H']))

    checkTrue(nrow(myFlowFrame[r1 == TRUE]) == 5206)
    checkTrue(nrow(myFlowFrame[r1 == FALSE]) == 14794)
    checkTrue(nrow(myFlowFrame[r2 == TRUE]) == 13752)
    checkTrue(nrow(myFlowFrame[r2 == FALSE]) == 6248)
}
