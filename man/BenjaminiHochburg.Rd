\name{BenjaminiHochburg}
\alias{BenjaminiHochburg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Finds significant features using false discovery rate (FDR).}
\description{
This function finds significant features using multiple test method proposed by Benjamini 
and Hochberg(1995). A threshold of FDR should be specified. 
}
\usage{
BenjaminiHochburg(p, q = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{A vector of p-values.}
  \item{q}{Specified threshold level. The default value is 0.05.}
}


\value{
Order of the most significant features corresponding to the selected threshold.
}
\references{ Benjamini, Y. and Hochberg, Y. (1995) Controlling the false discovery rate: a practical
and powerful approach to multiple testing. J. R. Statist. Soc. B, 57, 289-300. }

\author{ A. I. McLeod and M. S. Islam}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}