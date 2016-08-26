test_get_results_for_dyes_na <- function() {
    dyes <- c('a', 'b', 'c')
    detectors <- c('e', NA, 'e')
    results <- data.frame(x=1)
    checkTrue(nrow(get_results_for_dyes(dyes, detectors, results)) == 1)
}