\name{MeansBreakPoints}
\alias{MeansBreakPoints}

\title{ Finds the region means for given break points.}
\description{
 This function finds the region means when break points are give. 
This is one of the utility functions. 
}
\usage{
MeansBreakPoints(x, bp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Vector of observations }
  \item{bp}{Break points }
}

\value{
  Region means.
}

\references{ 
Islam, M. S. (2008). Periodicity, Change Detection and Prediction in Microarrays.
Ph.D. Thesis, The University of Western Ontario.
}

\author{ A. I. McLeod and M. S. Islam }


\examples{
set.seed(123)
x<-c(rnorm(20), rnorm(20, 1), rnorm(20))
MeansBreakPoints(x, c(20,30,50))
}

\keyword{nonparametric}
