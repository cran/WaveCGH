\name{PvalBootOneModel}
\alias{PvalBootOneModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Finds the p-values for testing means}
\description{
This function calculates p-values for testing means when we consider the
residuals obtained from all chromosomes to produce a single series. The break
means are obtained using wavelet transformation in each individual chromosomes. 
Then the p-values for testing mean are obtained using bootstrap procedure.  
}
\usage{
PvalBootOneModel(z, Chromosome, ModelSelection = FALSE, nsim = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ Array-CGH data. }
  \item{Chromosome}{ Sequence of chromosomes, eg. 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 and so on. }
  \item{ModelSelection}{ Whether to estimate order (p) of an AR(p) process.}
  \item{nsim}{Number of bootstrap simulations. }
}

\value{
Vector of p-values of length equal to that of z.
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
PvalBootOneModel(arrayCGH[(ch==1|ch==2),1], arrayCGH[(ch==1|ch==2),2])
}

\keyword{nonparametric}
