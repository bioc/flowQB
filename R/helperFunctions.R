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

find_peak <- function(data, width = .5, fraction = .1)
{
    x <- sort(data)
    N <- length(data)
    M <- ceiling(N * fraction)
    M2 <- floor((M + 1) / 2)
    for (first in 1:(N - M - 1))
        if (x[first] < x[first + 1])
            break
    
    for (last in N:(first + M))
        if (x[last - 1] < x[last])
            break
    
    i <- first
    dx <- x[last] - x[first]
    lo <- first
    hi <- last - M
    
    # Find the densest fraction (default 10%) of cells
    for (j in first:(last - M + 1))
    {
        d <- x[j + M - 1] - x[j]
        d
        if (d < dx)
        {
            i <- j
            dx <- d
        }
    }
    
    # Find the lower bound (default width at half maximum)
    for (j in i:first)
    {
        d <- x[j + M - 1] - x[j]
        if (d > dx / width)
        {
            lo <- j
            break
        }
    }
    # Find the upper bound (default width at half maximum)
    for (j in i:(last - M + 1))
    {
        d <- x[j + M - 1] - x[j]
        if (d > dx / width)
        {
            hi <- j
            break
        }
    }
    
    list(lo = x[lo + M2], hi = x[hi + M2])
}
