\name{etasFLP-package}
\alias{etasFLP-package}
\alias{etasFLP}
\docType{package}
\title{Mixed FLP and ML Estimation of ETAS Space-Time Point Processes
}
\description{
Estimation of the components of an ETAS model for earthquake description. 
 Non-parametric background seismicity can be estimated through FLP (Forward Likelihood Predictive), while parametric components are estimated through maximum likelihood. 
 The two estimation steps are alternated until convergence is obtained. 
 For each event the probability of being a background event is estimated and used as a weight for declustering steps. Many options to control the estimation process are present, together with some diagnostic tools. Some descriptive functions for earthquakes catalogs are  present; also \code{plot}, \code{print}, \code{summary}, \code{profile} methods are defined for main output (objects of class \code{etasclass})}

\details{
\tabular{ll}{
Package: \tab etasFLP\cr
Type: \tab Package\cr
Version: \tab 1.3.1\cr
Date: \tab 2016-12-07\cr
License: \tab GPL (>=2) \cr
Depends: \tab R (>= 2.14.0), mapdata,rgl \cr
Suggests: \tab MASS \cr
}
\code{etasclass} is the main function of the package \code{etasFLP}.
}
\author{
Marcello Chiodi and Giada Adelfio

Maintainer: Marcello Chiodi<marcello.chiodi@unipa.it>
}

\note{The package is intended for the estimation of the ETAS model for seismicity description (introduced by Ogata (1988), see reference),
but theoretically it can be used for other fields of application.}

\references{

Adelfio, G. and Chiodi, M. (2013) Mixed estimation technique in semi-parametric space-time point processes for earthquake description. 
\emph{Proceedings of the 28th International Workshop on Statistical Modelling 8-13 July, 2013, Palermo} (Muggeo VMR, Capursi V, Boscaino G, Lovison G, editors). Vol. \bold{1}. pp.65-70.

Adelfio, G. and Chiodi, M. (2015) Alternated estimation in semi-parametric space-time branching-type point processes with application to seismic catalogs. 
\emph{Stochastic Environmental Research and Risk Assessment}   DOI: 10.1007/s00477-014-0873-8 \bold{29}(2), pp. 443-450.


Chiodi, M. and Adelfio, G., (2011) Forward Likelihood-based predictive approach for space-time processes. \emph{Environmetrics}, vol. \bold{22} (6), pp. 749-757.

Console, R.,  Jackson,   D. D. and  Kagan, Y. Y.
Using the ETAS model for Catalog Declustering and Seismic
Background Assessment. \emph{Pure Applied Geophysics}
\bold{167}, 819--830 (2010).



Ogata, Y.   Statistical models for earthquake
occurrences and residual analysis for point processes.
\emph{Journal of the American Statistical Association},
\bold{83},  9--27 (1988).



Veen, A. and Schoenberg, F.P.  Estimation of space-time
branching process models in seismology using an EM-type algorithm.
\emph{Journal of the American Statistical Association}, \bold{103}(482), 614--624 (2008).


Zhuang, J., Ogata, Y.  and Vere-Jones, D.
Stochastic declustering of space-time earthquake occurrences.
\emph{Journal of the American Statistical Association},
\bold{97},  369--379 (2002).



}

\keyword{package}
\keyword{earthquakes}
\keyword{ETAS }
\keyword{PointProcess}
\keyword{MLE} 
\keyword{FLP} 
