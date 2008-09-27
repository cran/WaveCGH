\name{WaveCGH-package}
\alias{WaveCGH-package}
\alias{WaveCGH}
\docType{package}
\title{
Finding chromosomal gain/loss regions using wavelet method.
}
\description{
It finds and plots DNA copy number  gain/loss regions using wavelet method. Double implementation of MODWT confirms the presence of break points in the series. Bootstrapping procedure is used to call different regions to be normal, gain or loss regions. 
}
\details{
\tabular{ll}{
Package: \tab WaveCGH\cr
Type: \tab Package\cr
Version: \tab 2.01\cr
Date: \tab 2008-09-27\cr
License: \tab GPL (>= 2)\cr
}

}
\author{
M. S. Islam and A. I. McLeod
Maintainer: <aimcleod@uwo.ca>
}

\references{ 
Islam, M. S. (2008). Periodicity, Change Detection and Prediction in Microarrays.
Ph.D. Thesis, The University of Western Ontario.
}

\keyword{ package }

\examples{
data(arrayCGH)
# Finding break means for each chromosome
BreakMeansChr(arrayCGH[,1], arrayCGH[,2])
# Finding pvalues for testing means:
PvalBootMeanChr(arrayCGH[,1], arrayCGH[,2])
}
