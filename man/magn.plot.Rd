\name{magn.plot}
\alias{magn.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transformed plot of the magnitudes distribution of an earthquakes catalog}
\description{
Plots the logarithm of the cumulative frequency of eccedence vs. magnitude in an earthquake catalog.
}
\usage{
magn.plot(catalog,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{catalog}{should be a \code{eqcat} object, or at least must contain a column with name \code{magn1}.}
  \item{...}{other arguments to be passed to \code{plot()}}
}
\details{
For each magnitude \eqn{m_j}{mj}, if \eqn{N_j}{Nj} is the number of values of \code{magn1} greater than \eqn{m_j}{mj},
the values of  \eqn{\log(N_j)}{log(Nj)} vs. \eqn{m_j}{mj} are plotted.

According to the Gutenberg-Richter law, this plot should be linear. If there is a linear behaviour only for values greater 
than a given \eqn{m_0}{m0}, then  \eqn{m_0}{m0} is probably the magnitude threshold of the catalog.
}

\value{A new plot is printed (see details).}
%\references{}
\author{Marcello Chiodi.}

\examples{
\dontrun{
data(italycatalog)
magn.plot(italycatalog) 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{magnitude}
\keyword{Gutenberg-Richter}% __ONLY ONE__ keyword per line
