\name{PvalBootMeanChr}
\alias{PvalBootMeanChr}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Finds p-values for testing means in individual chromosomes. }
\description{
 This function applies maximum overlap discrete wavelet transform (MODWT) to find
break means inside each of the chromosomes. Then applies bootstrap method to find 
p-values for testing the means.  
}

\usage{
PvalBootMeanChr(x, Chromosome, ModelSelection = FALSE, nsim = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Array-CGH data.}
  \item{Chromosome}{ Sequence of chromosomes, eg. 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 and so on. }
  \item{ModelSelection}{ Whether to estimate order (p) of an AR(p) process.}
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
PvalBootMeanChr(arrayCGH[(ch==1|ch==2),1], arrayCGH[(ch==1|ch==2),2])
}

\keyword{nonparametric}
