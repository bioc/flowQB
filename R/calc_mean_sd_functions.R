###############################################################################
## Copyright (c) 2016
## Josef Spidlen, Faysal El Khettabi, Wayne Moore, David Parks, Ryan Brinkman
##
## License
## The software is distributed under the terms of the
## Artistic License 2.0
## http://www.r-project.org/Licenses/Artistic-2.0
## 
## Disclaimer
## This software and documentation come with no warranties of any kind.
## This software is provided "as is" and any express or implied
## warranties, including, but not limited to, the implied warranties of
## merchantability and fitness for a particular purpose are disclaimed.
## In no event shall the copyright holder be liable for any direct,
## indirect, incidental, special, exemplary, or consequential damages
## (including but not limited to, procurement of substitute goods or
## services; loss of use, data or profits; or business interruption)
## however caused and on any theory of liability, whether in contract,
## strict liability, or tort arising in any way out of the use of this
## software.
###############################################################################

calc_mean_sd_197 <- function(fcs_file_path, scatter_channels, ignore_channels)
{
    if (!file.exists(fcs_file_path)) return()
    result <- data.frame(
        row.names=c("total", "scatter gated", "peak gated", "mean", "sd"))
    fcs <- read.FCS(fcs_file_path)

    ## The fitted_ellipse_gate does that already, so the only effect this has
    ## is on the result total number of cells, which are eiher all cells, or
    ## just those cells with positive scatter values. 
    ## All cells is better I think.
    # for (scatter in scatter_channels)
    # {
    #     data <- exprs(fcs[,scatter])
    #     fcs <- fcs[data>0]
    # }
    scatter.gated <- fitted_ellipse_gate(fcs, scatter_channels, 2)
    fluorescences <- pick_parameters(fcs, ignore_channels)
    
    for (fl in fluorescences)
    {
        peak.gated <- scatter.gated[peak_gate(scatter.gated, fl, 4)]
        out <- getOutliers(as.vector(exprs(peak.gated[, fl])), 
            distribution="normal")
        result[[fl]] <- c(nrow(fcs), nrow(scatter.gated), 
            nrow(peak.gated), out$mu, out$sigma)
    }
    result
}

calc_mean_sd_duke <- function(fcs_file_path, scatter_channels, ignore_channels)
{
    ## The original calculation for duke and 197 beads is the same except that
    ## the 197 calculation putputs "total", "scatter gated", "peak gated", 
    ## "mean", and "sd" while the duke calculation outputs only "total", 
    ## "scatter gated", "mean", and "sd", and duke's "scatter gated" is actually
    ## "peak gated" in 197, which is more appropriate since it is created by a
    ## peak gate applied on the scatter gate. So we will drop the original
    ## duke calculation and use the 197 calculation instead.
    calc_mean_sd_197(fcs_file_path, scatter_channels, ignore_channels)

    # if (!file.exists(fcs_file_path)) return()
    # result <- data.frame(row.names=c("total", "scatter gated", "mean", "sd"))
    # fcs <- read.FCS(fcs_file_path)
    # scatter.gated <- fitted_ellipse_gate(fcs, scatter_channels, 2)
    # 
    # fluorescences <- pick_parameters(fcs, ignore_channels)
    # for (fl in fluorescences)
    # {
    #     peak.gated <- scatter.gated[peak_gate(scatter.gated, fl, 4)]
    #     out <- getOutliers(as.vector(exprs(peak.gated[, fl])), 
    #         distribution="normal")
    #     result[[fl]] <- c(nrow(fcs), nrow(peak.gated), out$mu, out$sigma)
    # }
    # result
}

calc_mean_sd_capture <- function(fcs_file_path, scatter_channels, detector, dye)
{
    if (!file.exists(fcs_file_path)) return()
    result <- data.frame(row.names=c(
        "total", "scatter gated", "high peak", "low peak", "stained mean",
        "stained sd", "unstained mean", "unstained sd"))
    fcs <- read.FCS(fcs_file_path)
    scatter.gated <- fitted_ellipse_gate(fcs, scatter_channels, 2)
    split <- split_in_two(scatter.gated, detector)
    capture.beads <- scatter.gated[split]
    capture.beads <- capture.beads[peak_gate(capture.beads, detector, 4)]
    s_out <- getOutliers(as.vector(exprs(capture.beads[, detector])),
        distribution="normal")
    blank.beads <- scatter.gated[!split]
    blank.beads <- blank.beads[peak_gate(blank.beads, detector, 4)]
    u_out <- getOutliers(as.vector(exprs(blank.beads[, detector])),
        distribution="normal")

    result[[dye]] <- c(nrow(fcs), nrow(scatter.gated), nrow(capture.beads),
        nrow(blank.beads), s_out$mu, s_out$sigma, u_out$mu, u_out$sigma);
    result
}

calc_mean_sd_capture_all <- function(fcs_file_path_list, scatter_channels_list,
    detector_list, dye_list)
{
    res <- apply(cbind
        (fcs_file_path_list, scatter_channels_list, detector_list, dye_list),
        1,
        function(x) calc_mean_sd_capture(x[[1]], x[[2]], x[[3]], x[[4]])
    )
    ret <- data.frame(matrix(unlist(res), nrow=nrow(res[[1]]), byrow=FALSE),
        row.names=rownames(res[[1]]))
    colnames(ret) <- dye_list
    ret
}

calc_mean_sd_background <- function(fcs_file_path, ignore_channels)
{
    if (!file.exists(fcs_file_path)) return()
    result <- data.frame(row.names=c("total", "mean", "sd"))
    fcs <- read.FCS(fcs_file_path)
    fluorescences <- pick_parameters(fcs, ignore_channels)

    for (fl in fluorescences)
    {
        out <- getOutliers(as.vector(exprs(fcs[, fl])), distribution="normal")
        result[[fl]] <- c(nrow(fcs), out$mu, out$sigma)
    }
    result
}
