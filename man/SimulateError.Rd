\name{SimulateError}
\alias{SimulateError}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulating AR(p) error series. }
\description{
This function finds an AR(P) error series.   
}
\usage{
SimulateError(phi, SigE, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{ Vector of AR(p) coefficients. }
  \item{SigE}{Standard deviation of errors in an AR(p) series.}
  \item{n}{ Length of series. }
}

\value{
Simulated series with AR(p) process.
}

\author{ A. I. McLeod and M. S. Islam }


\keyword{nonparametric}
