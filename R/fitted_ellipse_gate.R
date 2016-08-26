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

setGeneric(
    "fitted_ellipse_gate",
    def=function(object, channels, R=1) standardGeneric("fitted_ellipse_gate"),
    useAsDefault=function(object, channels, R=1)
    {
        stop(paste(
            "The fitted_ellipse_gate method is not supported on object type:",
            class(object)))
    }
)

## TODO: originally called scatter_gate
## This is applicable to any channels provided as arguments
## to this function, so it has been renamed, but the rest of the code still
## needs to be refactored accordingly
##
## Fit an ellipse (or ellipsoid) gate on the most dense region 
## of the selected channels of a flowFrame object.
setMethod(
    "fitted_ellipse_gate",
    signature=signature(object="flowFrame"),
    definition=function(object, channels, R=1)
    {
        ## Remove events with negative values in the selected channels 
        ## so that the log transform later doesn't fail
        for (channel in channels)
        {
            values <- exprs(object[,channel])
            object <- object[values>0]
        }
        
        ## Gate by an ellipse (ellipseoid) on the log transformed
        ## values of selected channels
        R2 <- vector(mode="double", nrow(object))
        for (channel in channels)
        {
            values <- log(exprs(object[,channel]))
            fwhm <- find_peak(values)
            center <- (fwhm$hi + fwhm$lo)/2
            radius <- (fwhm$hi - fwhm$lo)/2
            R2 = R2 + ((values - center)/radius)^2
        }
        
        ## Final gate can be stretched/shrunk by specifying R different from 1
        object <- object[R2 <= R^2]
        ## This does not fix the TOT keyword, so we will fix that manually
        description(object)$`$TOT` <- as.character(nrow(object))
        object
    }
)

