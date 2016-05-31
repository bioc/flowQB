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

calc_mean_sd_duke <- function(fcsFilePath, scatter_channels, ignore_channels)
{
    if (!file.exists(fcsFilePath)) return()
    result <- data.frame(row.names=c("total", "scatter gated", "mean", "sd"))
    fcs <- read.FCS(fcsFilePath)
    scatter.gated <- fitted_ellipse_gate(fcs, scatter_channels, 2)

    fluorescences <- pick_parameters(fcs, ignore_channels)
    for (fl in fluorescences)
    {
        peak.gated <- scatter.gated[peak_gate(scatter.gated, fl, 4)]
        out <- getOutliers(as.vector(exprs(peak.gated[, fl])), 
            distribution="normal")
        result[[fl]] <- c(nrow(fcs), nrow(peak.gated), out$mu, out$sigma)
    }
    result
    # dye.result <- get_results_for_dyes(xlsx, result)
}
