\name{plot.profile.etasclass}
\alias{plot.profile.etasclass}
%\alias{plot.etasclass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot method for profile.etasclass objects (profile likelihood of ETAS model)}
\description{\code{plot} method for \code{profile.etasclass} objects (profile likelihood of ETAS model). Plots a smooth interpolation of the profile likelihood of a parameter of an ETAS model, as output from \code{profile.etasclass}.
}
\usage{
\method{plot}{profile.etasclass}(x,prob=c(0.90,0.95,0.99),...)
}
\arguments{
  \item{x}{An object of the class  \code{profile.etasclass}.}
  \item{prob}{A vector of coverage probability for the asymptotic confidence interval computed using  -2log(LR). Default value \code{prob}=c(0.90,0.95,0.99).}
%%     ~~Describe \code{profile} here~~

  \item{...}{other arguments.}
}
\details{Plots a spline interpolation of the profile likelihood for a parameter of the ETAS model for earthquake seismicity, computed with \code{profile.etasclass}; 

the order of parametrs is: (\code{mu}, \code{k0}, \code{c}, \code{p}, \code{a}, \code{gamma}, \code{d}, \code{q}).

A \code{plot} method is defined for \code{profile.etasclass} objects. A number of grid points \code{nprofile}  of 7 (the default) usually is enough to have a good interpolation of the profile likelihood.
}

\value{Plots a profile likelihood (in the scale -2log(LR)), and plots horizontal lines corresponding to the percentiles of a 1df chi-square variable of levels \code{prob}; the approximate confidence intervals corresponding to the levels \code{prob} are printed. Returns a list:
\item{spline.profile}{The spline interpolation of the profile likelihood.}
\item{conf}{The approximate confidence intervals corresponding to the levels \code{prob}.}
\item{prob}{The \code{prob} values used.}
}
\note{A odd number of grid points \code{nprofile}  is adviced, so that the central point is the unconstrained ML estimate for the profiled parameter, and the interpolation of the profile likelihood will have a better quality.}

\author{
Marcello Chiodi, Giada Adelfio}

\seealso{ \code{\link{eqcat}}, \code{\link{etasclass}}, \code{\link{profile.etasclass}}}



\examples{
\dontrun{
## see example in profile.etasclass
}}
\keyword{ETAS}
\keyword{earthquake}
\keyword{profile}
\keyword{likelihood}

