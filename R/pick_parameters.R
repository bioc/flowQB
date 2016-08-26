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
    "pick_parameters",
    def=function(object, ignore) standardGeneric("pick_parameters"),
    useAsDefault=function(object, ignore)
    {
        stop(paste(
            "The pick_parameters method is not supported on object type:",
            class(object)))
    }
)

## TODO: originally called get_fluorescences
## This is applicable to picking any channels from a flowFrame, so it has been
## renamed, but the rest of the code still## needs to be refactored accordingly
setMethod(
    "pick_parameters",
    signature=signature(object="flowFrame"),
    definition=function(object, ignore)
    {
        ## Bellow is the original code, which I am replacing with a setdiff
        # pars <- colnames(object);
        # picked <- vector(mode='logical', length=length(pars))
        # for (j in 1:length(pars))
        # {
        #     picked[j] <- TRUE
        #     for (k in 1:length(ignore)) {
        #         if (pars[j] == ignore[k]) picked[j] <- FALSE
        #     }
        # }
        # pars[picked]
        
        ## New code:
        setdiff(colnames(object),ignore)
    }
)
