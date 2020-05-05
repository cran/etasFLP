\name{simpson.coeff}
\alias{simpson.coeff}
\alias{simpson.kD}
\title{
Computes Simpson integration rule coefficients}
\description{
Computes Simpson integration rule coefficients.}
\usage{
simpson.coeff(n)
simpson.kD(n,k=2)
}
\arguments{
  \item{n}{
number of points of the simpson formula a single dimension}
  \item{k}{
number of dimensions}
}

\details{
\code{simpson.coeff} computes the coefficients of the standard Simpson rule (for unit spaced points), according to the sequence (1+4+2+4+...+2+4+1)/3 for each dimension. \code{simpson.kD} expand the formula over a grid of \eqn{n^k} points in \code{k} dimensions.}

\value{
a vector of \code{n} coefficients (for \code{simpson.coeff}), a k-dimensions array with a total of \eqn{n^k} elements for \code{simpson.kD}.}

\author{
Marcello Chiodi
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{simpson}
\keyword{quadrature}