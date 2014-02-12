\name{b.guten}
\alias{b.guten}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates the parameter of the Gutenberg-Richter law.}
\description{
Estimates the parameter of the Gutenberg-Richter law for the magnitude distribution   of earthquakes, given a threshold magnitude; it
uses moment estimator on transformed data.}
\usage{
b.guten(magn, m0=min(magn))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{magn}{
a vector of magnitudes  coming from an earthquake catalog.
}
  \item{m0}{A threshold value. Only values of \code{magn} not less than \code{m0} will be used.}
}

\details{
Maximum likelihodd estimation for the Gutenberg-Richter Law:

\deqn{\log_{10} N(>m) = a - b \ M  }{log10 N(>m) = a - b m  }

where \eqn{N(>m)} is the number of events exceeding a magnitude \eqn{m} and \eqn{a,b} are two parameters: \eqn{a}  is related to the total seismicity rate of the region while \eqn{b}, to be estimated, should be usually near 1.

Catalog is assumed to be complete (in a certain space-time region) at least for a magnitude \code{m0},
that is, every earthquake of magnitude at least \code{m0} in that space-time region, is certainly present in the catalog.}
\value{
  \item{b}{estimate of the parameter  \eqn{b} of the Gutenberg-Richter Law.}
  \item{se}{estimate of the standard error of  the estimate \code{b}.}
}
\note{the plot produced by \code{magn.plot} can be used to have an idea, for a given catalog, of the magnitude threshold value.}
\seealso{ \code{\link{magn.plot}}}

\references{
Gutenberg, B. and Richter, C. F. (1944). Frequency of earthquakes in California.
\emph{Bulletin of the Seismological Society of America}, 34, 185-188.
}
\author{
Marcello Chiodi
}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{}
\examples{

data(italycatalog)
b.guten(italycatalog$magn1)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{magnitude}
\keyword{Gutenberg}
\keyword{Richter}
\keyword{earthquake}% __ONLY ONE__ keyword per line
