\name{BreakMeansChr}
\alias{BreakMeansChr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Finds break means in individual chromosomes. }
\description{
 This function applies maximum overlap discrete wavelet transform (MODWT) to find
break means inside each of the chromosomes. Double implementation of MODWT confirms 
the presence of break points in the series.    
}
\usage{
BreakMeansChr(x, Chromosome)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ CGH array.}
  \item{Chromosome}{Sequence of chromosome number,  eg. 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 and so on.}
}

\value{
Vector of break means of length equal to that of original series is returned.  
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
BreakMeansChr(arrayCGH[(ch==1|ch==2),1], arrayCGH[(ch==1|ch==2),2])
}

\keyword{nonparametric}
