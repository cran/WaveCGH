\name{PlotAllOneModel}
\alias{PlotAllOneModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plots all gain and loss regions in one one plot. }
\description{
This function finds break means for individual chromosome and then finds p-values by considering the residuals
from all chromosomes as single series. Thereafter, it plots all the gain and loss regions in one plot.  
}
\usage{
PlotAllOneModel(x, Chromosome, ModelSelection = FALSE, nsim = 200, threshold = 0.05, qvalue = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Array-CGH data.}

  \item{Chromosome}{ Sequence of chromosomes, eg. 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 and so on. }
  \item{ModelSelection}{ Whether to estimate order (p) of an AR(p) process.}
  \item{nsim}{Number of bootstrap simulations. }
  
\item{threshold}{ FDR threshold for multiple comparison.}
  \item{qvalue}{Whether to use Storey's q-value approach or not. }
  \item{\dots}{ Other arguments plot command. }
}

\value{
Plot of chromosome loss/gain regions.
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
PlotAllOneModel(arrayCGH[(ch==1|ch==2),1], arrayCGH[(ch==1|ch==2),2])
}

\keyword{nonparametric}
