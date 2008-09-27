\name{ACF.plot}
\alias{ACF.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot ACF for the fitted residuals }
\description{
 This function estimates the break means for each of the chromosomes, then finds
the residuals by subtracting the means from original values. It plots ACF for such
residuals.
}
\usage{
ACF.plot(x, cl, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Value of CGH array. }
  \item{cl}{ A vector referring to the sequence of chromosome numbers. It has same 
length as that of x.}
  \item{\dots}{ Other arguments that acf takes. }
}

\value{
No value is returned, only ACF plot. 
}

\author{ A. I. McLeod and M. S. Islam}

\seealso{ \code{\link{acf}}}

\examples{
set.seed(123)
x<-c(rnorm(20), rnorm(20, 1), rnorm(20))
cl<-rep(1, length(x)) # one chromosome
ACF.plot(x, cl)
}

\keyword{nonparametric}
