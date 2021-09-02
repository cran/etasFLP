\name{profile.etasclass}
\alias{profile.etasclass}
%\alias{plot.etasclass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{profile method for etasclass objects (ETAS model) (To be checked)}
\description{ profile method for etasclass objects (ETAS model).
}
\usage{
\method{profile}{etasclass}(fitted,iprofile		=1,
				 nprofile	=7,
				 kprofile	=3,
				 profile.approx	=FALSE,...)
}
\arguments{
  \item{fitted}{An object of the class  \code{etasclass}}.
  \item{iprofile}{An integer in the range 1-7. Profile likelihood will be computed with respect to the parameter of index \code{iprofile}. The order of parametrs is: \code{mu}, \code{k0}, \code{c}, \code{p}, \code{a}, \code{gamma}, \code{d}, \code{q}.
In this version 2.0.0 \code{profile} only for etas parameters, not for covariates parameters. It will added future versions  				 
}
  \item{nprofile}{Number of values of \code{params[iprofile]} for which profile likelihood must be computed.  Default value= 7.
%%     ~~Describe \code{nprofile} here~~
}
  \item{kprofile}{Maximum absolute standardized value for \code{params[iprofile]}. Profile likelihood will be computed in the standardized range [-\code{kprofile}, \code{kprofile}].  Default value= 3.
%%     ~~Describe \code{kprofile} here~~
}
  \item{profile.approx}{if \code{TRUE}, then a conditional-likelihood approach is used as a first value for  each maximization step in profile likelihood computation.  Default value= \code{FALSE}.
%%     ~~Describe \code{profile.approx} here~~
}
  \item{...}{other arguments.}  
}
\details{Profile likelihood for the \code{iprofile}-th parameter of the ETAS model for earthquake seismicity, estimated with \code{etasclass}; the order of parameters is: \code{mu},\code{k0},\code{c},\code{p},\code{gamma},\code{d},\code{q} and \code{betacov}.

A \code{plot} method is defined for \code{profile.etasclass} objects. A number of grid points \code{nprofile}  of 7 (the default) usually is enough to have a good interpolation of the profile likelihood. The profile is computed using the final estimation of the background seismicity used to obtain the object \code{etas} of class \code{etasclass} and regardless to the method used. The computing time (for each of the \code{nprofile} values) is generally less than a single execution of \code{etasclass} without clustering, because only ML estimation is performed. Parameters not estimated in \code{etas} (with \code{params.ind[i]=FALSE}) will remain fixed do the value \code{params.fix[i]}.

To obtain profiles for different parameters, run \code{profile.etasclass} with different values of \code{iprofile}.
}

\value{Returns a list:
\item{params.vec}{vector of values of the parameter \code{iprofile} used to evaluate the profile likelihood.}
\item{logl.vec}{vector of likelihoods corresponding to the values of \code{params.vec}}
\code{plot} method is defined to represent  profile likelihood (in  scale -2log(LR)), using a spline interpolation through  grid points, with superimposition of approximate confidence intervals.
}
\note{A odd number of grid points \code{nprofile}  is adviced, so that the central point is the unconstrained ML estimate for the profiled parameter, and the interpolation of the profile likelihood will have a better quality.}

\author{
Marcello Chiodi, Giada Adelfio}

\seealso{ \code{\link{eqcat}}, \code{\link{etasclass}}, \code{\link{plot.profile.etasclass}}}



\examples{
\dontrun{##
data("italycatalog")
# load a sample catalog of italian seismicity  

   etas.flp<-etasclass(italycatalog,  
   magn.threshold = 3.0,  magn.threshold.back = 3.5,
   k0 = 0.005, c = 0.005, p = 1.01, gamma = 0.6, q = 1.52, d = 1.1,
   params.ind = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
   formula1 = "time ~  magnitude- 1", declustering = TRUE, 
   thinning = FALSE, flp = TRUE, ndeclust = 15,
   onlytime = FALSE, is.backconstant = FALSE,
   description = "etas flp",sectoday = TRUE,
   usenlm = TRUE, epsmax = 0.001)
    # execution of etasclass for events with minimum magnitude of 3.0. 
    # The events with magnitude at least 3.5 are used to build a first approximation
    # for the background intensity function
    # (magn.threshold.back=3.5)
    
## compute profile likelihood for the first parameter (mu)
 system.time( prof.flp <- profile(etas.flp, nprofile = 7, iprofile = 1))
 plot(prof.flp)
#### output:
Asymptotic  confidence  intervals: 
  Coverage Lower Upper
1     0.90 0.335 0.376
2     0.95 0.334 0.378
3     0.99 0.329 0.385

}}
\keyword{ETAS}
\keyword{earthquake}
\keyword{profile}
\keyword{likelihood}

