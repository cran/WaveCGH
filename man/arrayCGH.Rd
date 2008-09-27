\name{arrayCGH}
\alias{arrayCGH}
\docType{data}
\title{ Array CGH data}
\description{
In CGH array, 2400 BAC clones were measured each with three replicates (Snijders
et al., 2001). Measurements for log base 2 intensity ratio are provided.   
}
\usage{data(arrayCGH)}
\format{
  A data frame with 1040 observations on the following 2 variables.
  \describe{
    \item{\code{value}}{a numeric vector}
    \item{\code{chromosome}}{a numeric vector}
  }
}

\source{
Snijders, A. M., Nowak, N., Segraves, R., Blackwood, S., Brown N., Conroy, J.,
Hamilton, G., Hindle, A. K, Huey, B., Kimura, K., Law, S., Myambo, K., Palmer,
J., Ylstra, B., Yue, J. P., Gray, J. W., Jain, A. N., Pinkel, D. and Albertson,
D. (2001)Assembly of microarrays for genome-wide measurement of DNA copy
number. Nature Genetics 29, 263 - 264. 
}

\examples{
data(arrayCGH)
names(arrayCGH)
}
\keyword{datasets}
