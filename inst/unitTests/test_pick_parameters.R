test_pick_parameters <- function() {
    library('flowCore')
    fcsFilePath <- system.file("extdata", "SSFF_LSRII", "Other_Tests",
        "933745.fcs", package="flowQBData")
    myFlowFrame <- read.FCS(fcsFilePath)
    ignore <- c("Time", "FSC-H", "FSC-A", "FSC-W", "SSC-H", "SSC-A", "SSC-W",
        "foo", "doo")
    fluorescences <- pick_parameters(myFlowFrame, ignore)
    
    checkTrue("FITC-A" %in% fluorescences)
    checkTrue("FITC-H" %in% fluorescences)
    checkTrue("PerCP-Cy5-5-A" %in% fluorescences)
    checkTrue("PerCP-Cy5-5-H" %in% fluorescences)
    checkTrue("Pacific Blue-A" %in% fluorescences)
    checkTrue("Pacific Blue-H" %in% fluorescences)
    checkTrue("Aqua Amine-A" %in% fluorescences)
    checkTrue("Aqua Amine-H" %in% fluorescences)
    checkTrue("Pacific Orange-A" %in% fluorescences)
    checkTrue("Pacific Orange-H" %in% fluorescences)
    checkTrue("QDot 585-A" %in% fluorescences)
    checkTrue("QDot 585-H" %in% fluorescences)
    checkTrue("QDot 605-A" %in% fluorescences)
    checkTrue("QDot 605-H" %in% fluorescences)
    checkTrue("QDot 655-A" %in% fluorescences)
    checkTrue("QDot 655-H" %in% fluorescences)
    checkTrue("QDot 705-A" %in% fluorescences)
    checkTrue("QDot 705-H" %in% fluorescences)
    checkTrue("QDot 800-A" %in% fluorescences)
    checkTrue("QDot 800-H" %in% fluorescences)
    checkTrue("APC-A" %in% fluorescences)
    checkTrue("APC-H" %in% fluorescences)
    checkTrue("APC-Cy5-5-A" %in% fluorescences)
    checkTrue("APC-Cy5-5-H" %in% fluorescences)
    checkTrue("APC-Cy7-A" %in% fluorescences)
    checkTrue("APC-Cy7-H" %in% fluorescences)
    checkTrue("PE-A" %in% fluorescences)
    checkTrue("PE-H" %in% fluorescences)
    checkTrue("PE-Texas-Red-A" %in% fluorescences)
    checkTrue("PE-Texas-Red-H" %in% fluorescences)
    checkTrue("PE-Cy5-A" %in% fluorescences)
    checkTrue("PE-Cy5-H" %in% fluorescences)
    checkTrue("PE-Cy5-5-A" %in% fluorescences)
    checkTrue("PE-Cy5-5-H" %in% fluorescences)
    checkTrue("PE-Cy7-A" %in% fluorescences)
    checkTrue("PE-Cy7-H" %in% fluorescences)

    checkTrue(all(ignore %in% fluorescences == rep(FALSE, length(ignore))))
}

test_pick_parameters2 <- function() {
    expected <- "expected"
    tmp <- tryCatch(
        {
            pick_parameters("foo");
        },
        error = function(ex) {
            expected;
        }
    )
    checkTrue(tmp == expected)
}