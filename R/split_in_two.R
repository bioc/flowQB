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
    "split_in_two",
    def=function(object, ...) standardGeneric("split_in_two"),
    useAsDefault=function(object, ...)
    {
        stop(paste("The split_in_two method is not supported on object type:",
                   class(object)))
    }
)

setMethod(
    "split_in_two",
    signature=signature(object="matrix"),
    definition=function(object, ...)
    {
        x <- sort(object)
        N <- length(x)
        M <- ceiling(N/10)
        M2 <- floor((M + 1)/2)
        m <- i <- floor((N + 1)/2)
        dx <- x[m+M2] - x[m-M2];
        for (j in 1:ceiling(N/4))
        {
            d <- x[m+j+M2] - x[m+j-M2];
            if (d > dx)
            {
                dx <- d
                i <- m + j
            }
            d <- x[m-j+M2] - x[m-j-M2];
            if (d > dx)
            {
                dx <- d
                i <- m - j
            }
        }
        object >= x[i]
    }
)

setMethod(
    "split_in_two",
    signature=signature(object="flowFrame"),
    definition=function(object, channel, ...)
    {
        split_in_two(exprs(object[,channel]))
    }
)
