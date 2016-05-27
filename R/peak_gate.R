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
    "peak_gate",
    def=function(object, ...) standardGeneric("peak_gate"),
    useAsDefault=function(object, ...)
    {
        stop(paste("The peak_gate method is not supported on object type:",
                   class(object)))
    }
)

setMethod(
    "peak_gate",
    signature=signature(object="flowFrame"),
    definition=function(object, channel, R=1, ...)
    {
        peak_gate(exprs(object[,channel]), R=R)
    }
)

setMethod(
    "peak_gate",
    signature=signature(object="matrix"),
    definition=function(object, R=1, ...)
    {
        fwhm <- find_peak( object )
        center <- (fwhm$hi + fwhm$lo)/2
        radius <- (fwhm$hi - fwhm$lo)/2
        abs(object - center) <= R * radius
        
    }
)
