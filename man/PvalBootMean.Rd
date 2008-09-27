\name{PvalBootMean}
\alias{PvalBootMean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Finds p-values for testing means using bootstrap method.}
\description{
This function finds break means using MODWT. Then calculates p-values for testing means 
using bootstrap method.   
}
\usage{
PvalBootMean(x, ModelSelection = FALSE, nsim = 200)
}
\arguments{
  \item{x}{ Array-CGH data.}
  \item{ModelSelection} { Whether to estimate order (p) of an AR(p) process.}
  \item{nsim}{Number of bootstrap simulations. }
  
}

\value{
Vector of p-values of length equal to that of original series.
}

\references{ 
Islam, M. S. (2008). Periodicity, Change Detection and Prediction in Microarrays.
Ph.D. Thesis, The University of Western Ontario.
}

\author{ A. I. McLeod and M. S. Islam }


\examples{
data(arrayCGH)
ch<- arrayCGH[,2]
# We use first and second chromosomes
set.seed(123)
PvalBootMean(arrayCGH[(ch==1|ch==2),1])
}

\keyword{nonparametric}
