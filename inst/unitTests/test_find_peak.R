test_find_peak <- function() {
    set.seed(123456)
    peak_details <- find_peak(rnorm(1000, mean=5))
    checkTrue(length(peak_details) == 2)
    checkTrue(peak_details$lo > 4.03583)
    checkTrue(peak_details$lo < 4.03584)
    checkTrue(peak_details$hi > 6.08955)
    checkTrue(peak_details$hi < 6.08956)
}