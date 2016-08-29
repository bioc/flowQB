# test_split_in_two <- function() {
#     library('flowCore')
#     fcsFilePath <- system.file("extdata", "SSFF_LSRII", "Other_Tests",
#         "933745.fcs", package="flowQBData")
#     myFlowFrame <- read.FCS(fcsFilePath)
#     r1 <- split_in_two(myFlowFrame, 'FSC-H')
#     r2 <- split_in_two(exprs(myFlowFrame[,'SSC-H']))
# 
#     checkTrue(nrow(myFlowFrame[r1 == TRUE]) == 74998)
#     checkTrue(nrow(myFlowFrame[r1 == FALSE]) == 25002)
#     checkTrue(nrow(myFlowFrame[r2 == TRUE]) == 74870)
#     checkTrue(nrow(myFlowFrame[r2 == FALSE]) == 25130)
# }

test_split_in_two2 <- function() {
    expected <- "expected"
    tmp <- tryCatch(
        {
            split_in_two("foo", 'FSC-H');
        },
        error = function(ex) {
            expected;
        }
    )
    checkTrue(tmp == expected)
}