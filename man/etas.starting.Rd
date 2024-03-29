\name{etas.starting}
\alias{etas.starting}

\title{Guess starting values of ETAS parameters (beta-version). Only from package version 1.2.0}
\description{\code{etas.starting} is a simple function to give starting values of the 7 ETAS parameters for the function \code{etasclass}.

It gives only rough approximations, based on some assumptions, intended to give only the order of magnitude of each parameter (but should be better than nothing).
Returns a list with starting values. In the present version user can give manually the output of this function in the input of \code{etasclass}. Otherwise, the function is called by \code{etasclass} at first steps, to supply initial values to start estimation.
}
\usage{etas.starting(cat.orig,
magn.threshold=2.5,
p.start=1,
gamma.start=0.5,
q.start=2,
betacov.start=.7,
longlat.to.km=TRUE,
sectoday=FALSE,
onlytime=FALSE
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cat.orig}{
An  earthquake catalog, possibly an object of class \code{eqcat}, or however a \code{data.frame} with variables of names \code{time},  \code{lat},   \code{long},  \code{z},     \code{magn1}. No missing values are allowed.
}
\item{magn.threshold}{
Threshold magnitude (only events with a magnitude at least \code{magn.threshold} will be used).  Default value = 2.5.
}
\item{p.start}{Parameter 4  of the ETAS model; the exponent of the Omori law for temporal decay rate of aftershocks; see details. Default value = 1.0.
%%     ~~Describe \code{p} here~~
}
  \item{gamma.start}{Parameter 5 (\eqn{\gamma}) of the ETAS model; together with  \code{a} is related to the efficiency of an event of given magnitude 
  in generating aftershocks; see details. Default value = 0.5.
%%     ~~Describe \code{gamma} here~~
}
  \item{q.start}{Parameter 7 of the ETAS model; parameter related to the spatial influence of the
mainshock; see details. Default value = 2.
%%     ~~Describe \code{q} here~~
}
  \item{betacov.start}{coefficient of the covariate (as default the magnitude). Default value = 0.7.
%%     ~~Describe \code{q} here~~
}
  \item{sectoday}{ if \code{TRUE}, then \code{time} variable of \code{cat.orig} is converted from seconds to days.  Default value = \code{TRUE}.
}
  \item{longlat.to.km}{ if \code{TRUE}, then \code{long} and \code{lat} variables of \code{cat.orig} are treated as geographical coordinates and converted to kilometers.  Default value = \code{TRUE}.
}
\item{onlytime}{if \code{TRUE}  then a time process is fitted to data , regardless to space location (in this case \code{is.backconstant} is set to \code{TRUE} and \code{declustering}, \code{flp} are set to \code{FALSE}). Default value = \code{FALSE}.
}
}
\details{It is a beta-version of a very crude method to give 
 starting values for the seven parameters of an ETAS  (Epidemic type aftershock sequences) model
 for the description of the seismicity of a space-time region. 
These starting values can be used as input for the function \code{etasclass} 
\code{sectoday} and \code{ longlat.to.km} flags must the same that will be used in \code{etasclass}.

In this first attempt to give starting values for the ETAS model, many approximations are used 

It gives only rough approximation, based on some assumptions, intended to give only the order of magnitude of each parameter (but it should be better than nothing). It
returns a list with 7 starting values. With this beta-version user must give manually the output of this function in the input of \code{etasclass}. 

The values of \code{p.start},  \code{gamma.start} and  \code{q.start} must be however given by the user (we did not find anything reasonable). Default choices for \code{p} and \code{q} (\code{p.start=1}, \code{q.start=2}) are strongly reccomended.

\code{c} and \code{d} are estimated from the emprical distributions of time differences and space distances, respectively.
\code{mu} and \code{k0} are then estimated given the other starting values, solving the two ML equations, that is derivatives of the whole likelihood with respect to  \code{mu} and \code{k0} equated to zero. 
In the computation of the likelihood an approximation for the integral of the intensity function is used (quoted also in 
Schoenberg (2013)).
}

\value{
returns a list:
\item{mu.start}{guess value for \code{mu}}
\item{k0.start}{guess value for \code{k0}}
\item{c.start}{guess value for \code{c}}
\item{p.start}{guess value for \code{p} (the same as input value)}
\item{gamma.start}{guess value for \code{gamma} (the same as input value)}
\item{d.start}{guess value for \code{d}}
\item{q.start}{guess value for \code{q} (the same as input value)}
\item{longlat.to.km}{\code{longlat.to.km} (the same as input value)}
\item{sectoday}{\code{sectoday} (the same as input value)}
}

\note{
The optimization algorithm used in \code{etasclass}  depends on the choice of initial values. Some default guess choice is performed in the present beta-version of the function \code{etas.starting}. If convergence problem are experienced, a useful strategy can be to start  with an high magnitude threshold value \eqn{m_0}{m0} (that is, with a smaller catalog with bigger earthquakes), and then using this first output as starting  guess for a running with a lower magnitude threshold value \eqn{m_0}{m0}. 
In this trial executions avoid declustering (\code{declustering=FALSE}) or at least use a small value of \code{ndeclust}; small values of \code{iterlim} and \code{ntheta} can speed first executions.

Quicker executions are obtained using smaller values of \code{iterlim} and \code{ntheta} in the input.

Also a first execution with \code{is.backconstant = TRUE}, to fit a first approximation model with constant background, can be useful.

Some other useful information can be obtained estimating a pure time process, that can give a good guess at least for some  parameters, like  \eqn{\mu, \kappa_0, c,p}{mu, k0, a, c, p}.

Input times are expected in days, and so final intensities are expected number of events per day. If input values are in seconds, then set \code{sectoday=TRUE}
}

\seealso{ \code{\link{etasclass}}}

\references{
Schoenberg, F. P. (2013).Facilitated Estimation of ETAS. \emph{Bulletin of the Seismological Society of America}, Vol. 103, No. 1, pp. 601-605, February 2013, doi: 10.1785/0120120146
}

\author{
Marcello Chiodi, Giada Adelfio}

%% ~Make other sections like Warning with \section{Warning }{....} ~


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ETAS}
\keyword{earthquake}
\keyword{likelihood}
