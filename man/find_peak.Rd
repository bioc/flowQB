\name{find_peak}
\alias{find_peak}

\title{ 
    Find density peak in the provided vector of values
}

\description{
    Find a density peak in the provided vector of values and return
    the lower and upper bounds around that peak.
}

\usage{
    find_peak(data, width=0.5, fraction=0.1)
}

\arguments{
    \item{data}{
        Vector of numbers.
    }
    \item{width}{
        The returned lower and upper bounds for the peak are calculated
        based on how much larger is the difference between the values 
        at the "i"" and "i + faction size" index of sorted data is; i is
        itterated from the index of the most dense region to the left
        (for lower bound) and to the right (for upper bound), which increases
        the above mentioned difference as the region is becomming less
        dense. The itteration is stopped when the difference gets larger
        than the minimum difference times 1/width. This means, by default with
        width = 0.5, the difference can reach double the minimum difference.
    }
    \item{fraction}{
        How large is the proportion of values that we look at when distances
        as mentioned above. By default with fraction = 0.1, we look at 10\% 
        of values.
    }
}

\details{
    This functions finds a density peak in the provided vector of values and
    returns the lower and upper bounds around that peak. First, we sort the
    input data, let's call it x, and look at a fraction of values at a time.
    Say we have 1,000 values in our input vector and fraction is 0.1 (10\%),
    meaning we look at regions spanning 100 values. We will find the most 
    dense region by itterating the possible starting index from 1 to 899 and
    minimizing the difference between the sorted values at potential 
    starting index and corresponding ending index (100 values apart).
    Finding the minimum difference equals to finding the most dense region.
    Next we find the lower and upper bounds for the peak based on how much 
    larger is the difference between the values at the "i"" and 
    "i + faction size" index of sorted data is; i is itterated from the index
    of the most dense region to the left (for lower bound) and to the right
    (for upper bound), which increases the above mentioned difference as the
    region is becomming less dense. The itteration is stopped when the
    difference gets larger than the minimum difference times 1/width.
    This means, by default with width = 0.5, the difference can reach double
    the minimum difference.
}

\value{
    A list with 2 components, "lo" and "hi", identifying the lower and 
    the upper bounds for the peak.
}

\author{
    Wayne Moore, Faysal El Khettabi, Josef Spidlen
}

\seealso{
    \code{\link{fitted_ellipse_gate}}
}

\examples{
    my_peak_info <- find_peak(rnorm(1000, mean=5))
}

\keyword{functions}
\keyword{helper functions}
