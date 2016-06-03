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

fit_led <- function(fcs_file_path_list, ignore_channels, bounds, signal_type, 
    instrument_name)
{
    # TODO Rename peak.list to flowFrame_list
    peak.list <- list()
    # TODO Rename parameters to fluorescences
    parameters <- c()
    
    for (i in 1:length(fcs_file_path_list))
    {
        fcs <- read.FCS(fcs_file_path_list[[i]])
        peak.list[[i]] <- fcs
    }
    if (length(peak.list) > 0) 
        parameters <- pick_parameters(fcs, ignore_channels)
    
    results.list <- get_peak_statistics(peak.list, 
        basename(fcs_file_path_list), parameters, bounds, 
        signal_type == "Height", instrument_name == "BD Accuri")
    
    
    
}