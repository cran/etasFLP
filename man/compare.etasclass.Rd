\name{compare.etasclass}
\alias{compare.etasclass}

\title{
Compare two etasclass objects}
\description{
Compare the results of two \code{etasclass} executions through the comparison of some elements of two \code{etasclass} objects.}
\usage{
compare.etasclass(etas1,etas2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{etas1}{an \code{etasclass} object.}
  \item{etas2}{an \code{etasclass} object.}
}
\details{
Compare the results of two \code{etasclass} executions through the comparison of some elements of two \code{etasclass} objects. Fundamental elements (catalog, threshold magnitude) should be equal, but no check is made
}
\value{
\item{params}{Standardized comparison of the estimated parameters.}
\item{AIC}{Difference of AIC values.}
\item{weigths}{Comparison of the weights \code{rho.weights} of the two input objects.}
\item{rho.weigths}{Standardized comparison of the weights \code{rho.weights} of the two input objects.}
\item{cor.weights}{Correlation between the weights \code{rho.weights} of the two input objects}
\item{cor.trig}{Correlation between the triggered intensities of the two input objects}
\item{cor.back}{Correlation between the background intensities  of the two input objects}
}
\seealso{ \code{\link{etasclass}}}
\author{
Marcello Chiodi}

%\seealso{}
\keyword{ catalog }
\keyword{ etasclass }
\keyword{ earthquake }

