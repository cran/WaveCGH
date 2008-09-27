\name{BreakMeans}
\alias{BreakMeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Finds break means in a series using wavelet method. }
\description{
This function applies maximum overlap discrete wavelet transform (MODWT) to find
break means in a series. Double implementation of MODWT confirms the presence of
break points in the series. 
}
\usage{
BreakMeans(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{CGH array. }
}

\value{
Vector of break means corresponding to each point is returned.  
}


\references{ 
Islam, M. S. (2008). Periodicity, Change Detection and Prediction in Microarrays.
Ph.D. Thesis, The University of Western Ontario.
}

\author{ A. I. McLeod and M. S. Islam }


\examples{
set.seed(123)
x<-c(rnorm(20), rnorm(20, 1), rnorm(20))
BreakMeans(x)
}

\keyword{nonparametric}
