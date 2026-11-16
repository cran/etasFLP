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
Compare the results of two \code{etasclass} executions through the comparison of some elements of two \code{etasclass} objects estimated on the same catalog. If some standard erros are missing (or zero for constant parameters) they are not used in the standardized comparison of parameters (output object \code{diffstd.params}).
}
\value{
\item{diffstd.params}{Standardized comparison of the estimated parameters, included covariates.}
\item{AIC}{Difference of AIC values.}
\item{weights}{Comparison of the weights \code{rho.weights} of the two input objects.}
\item{weights.std}{Standardized comparison of the weights \code{rho.weights} of the two input objects.}
\item{cor.weights}{Correlation between the weights \code{rho.weights} of the two input objects}
\item{cor.trig}{Correlation between the triggered intensities of the two input objects}
\item{cor.back}{Correlation between the background intensities  of the two input objects}
}
\seealso{ \code{\link{etasclass}}}
\author{
Marcello Chiodi}

%\seealso{}
\examples{
data("italycatalog")
## a small sample of the catalog keeps the execution time short:
## please, use larger values of n to obtain reasonable datasets and estimates
sub <- italycatalog[1:150, ]
fit <- etasclass(sub, magn.threshold = 3.0, magn.threshold.back = 3.5,
    k0 = 0.005, c = 0.005, p = 1.01, gamma = 0.6, q = 1.52, d = 1.1,
    ndeclust = 1, sectoday = TRUE, epsmax = 0.001, trace = FALSE)

## a second object, obtained with a different number of iterations
fit2 <- update(fit, iterlim = 3)
compare.etasclass(fit, fit2)
}

\keyword{ catalog }
\keyword{ etasclass }
\keyword{ earthquake }

