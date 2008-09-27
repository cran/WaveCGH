\name{PlotLattice}
\alias{PlotLattice}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Plots all gain and loss regions in separate plots for chromosomes. }
\description{
This function finds break means for individual chromosome and then finds p-values for testing means. 
Thereafter, it plots all the gain and loss regions in separate plots.  
}

\usage{
PlotLattice(x, Chromosome, ModelSelection = FALSE, nsim = 200, threshold = 0.05, qvalue = FALSE, Panel = c(1, 1), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{x}{ Array-CGH data.}
  \item{Chromosome}{ Sequence of chromosomes, eg. 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 and so on. }
  \item{ModelSelection}{ Whether to estimate order (p) of an AR(p) process.}
  \item{nsim}{Number of bootstrap simulations. }
  
\item{threshold}{ FDR threshold for multiple comparison.}
  \item{qvalue}{Whether to use Storey's q-value approach or not. }
  \item{Panel}{The format of plotting the chromosomes.}
  \item{\dots}{ Other arguments plot command. }
}

\value{
Plot of gain/loss regions for each chromosome.
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
PlotLattice(arrayCGH[(ch==1|ch==2),1], arrayCGH[(ch==1|ch==2),2], Panel=c(1,2))
}

\keyword{nonparametric}
