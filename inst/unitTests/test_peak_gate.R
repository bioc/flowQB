test_peak_gate <- function() {
    library('flowCore')
    fcsFilePath <- system.file("extdata", "935289.fcs", package="flowQB")
    myFlowFrame <- read.FCS(fcsFilePath)

    r <- peak_gate(myFlowFrame, 'FSC-H')
    checkTrue(nrow(myFlowFrame[r == TRUE]) == 14133)
    checkTrue(nrow(myFlowFrame[r == FALSE]) == 5867)

    r2 <- peak_gate(myFlowFrame, 'FSC-H', 2)
    checkTrue(nrow(myFlowFrame[r2 == TRUE]) == 19034)
    checkTrue(nrow(myFlowFrame[r2 == FALSE]) == 966)
}
